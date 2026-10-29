import numpy as np
import matplotlib
matplotlib.use('Agg')
import lime
import pytest
import pyneb as pn
from pathlib import Path

from specsy.models.extinction import spectrum_reddening_correction, DEFAULT_RV, FIXED_RV_LAWS, EXT_LAWS, HBETA_WAVE
from specsy.io import SpecSyError


# Data for the tests (same baseline spectrum as test_spectrum.py)
baseline_folder = Path(__file__).parent / 'baseline'
file_address = baseline_folder/'SHOC579_MANGA38-35.txt'

redshift = 0.0475
norm_flux = 1e-17
res_power = 3000
e_bv_test = 0.25

wave_array, flux_array, err_array, pixel_mask = np.loadtxt(file_address, unpack=True)
pixel_mask = pixel_mask.astype(bool)

spec = lime.Spectrum(wave_array, flux_array, err_array, redshift=redshift, norm_flux=norm_flux,
                     pixel_mask=pixel_mask, res_power=res_power, id_label='SHOC579_Manga38-35')


def effective_rv(law, r_v=None):
    """R_V that PyNeb ends up using: the tabulated one by default and always for the fixed-R_V laws"""
    if (r_v is None) or (law in FIXED_RV_LAWS):
        return DEFAULT_RV[law]
    return r_v


def reference_klam(wave, law, r_v=None):
    """A_lambda/E(B-V) computed directly with PyNeb"""
    return np.squeeze(pn.RedCorr(E_BV=1.0, R_V=effective_rv(law, r_v), law=law).X(wave)).astype(float)


def reference_corr(wave, e_bv, law, r_v=None):
    return np.power(10, 0.4 * e_bv * reference_klam(wave, law, r_v))


def physical_flux(spectrum):
    return spectrum.flux.data * spectrum.norm_flux


def correction_ratio(spec_out, spec_in):
    """Output/input flux ratio, only meaningful (and returned) at the good pixels with non-zero flux"""
    flux_in, flux_out = physical_flux(spec_in), physical_flux(spec_out)
    idcs_good = ~spec_out.flux.mask & ~spec_in.flux.mask & np.isfinite(flux_in) & (flux_in != 0)
    return idcs_good, flux_out[idcs_good] / flux_in[idcs_good]


class TestSpectrumReddening:

    def test_input_spectrum_unchanged(self):

        wave_0, wave_rest_0 = spec.wave.data.copy(), spec.wave_rest.data.copy()
        flux_0, err_0 = spec.flux.data.copy(), spec.err_flux.data.copy()
        mask_0 = spec.flux.mask.copy()

        spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test)

        assert np.array_equal(wave_0, spec.wave.data)
        assert np.array_equal(wave_rest_0, spec.wave_rest.data)
        assert np.array_equal(flux_0, spec.flux.data, equal_nan=True)
        assert np.array_equal(err_0, spec.err_flux.data, equal_nan=True)
        assert np.array_equal(mask_0, spec.flux.mask)

        return

    def test_output_spectrum_properties(self):

        spec_out = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test)

        # New object of the same class
        assert isinstance(spec_out, lime.Spectrum)
        assert spec_out is not spec

        # Metadata
        assert spec_out.label == spec.label
        assert spec_out.redshift == redshift
        assert spec_out.norm_flux == norm_flux
        assert spec_out.res_power == res_power
        assert spec_out.units_wave == spec.units_wave
        assert spec_out.units_flux == spec.units_flux

        # The wavelength arrays do not change
        assert np.allclose(wave_array, spec_out.wave.data)
        assert np.allclose(wave_array / (1 + redshift), spec_out.wave_rest.data)

        # The measurements of the original spectrum are not carried over
        assert len(spec_out.frame) == 0
        assert spec_out.cont is None
        assert spec_out.cont_std is None

        return

    @pytest.mark.parametrize('rest_frame', [False, True])
    @pytest.mark.parametrize('law', EXT_LAWS)
    def test_flux_values(self, law, rest_frame):

        spec_out = spectrum_reddening_correction(spec, law=law, E_BV=e_bv_test, rest_frame=rest_frame)

        wave_law = wave_array / (1 + redshift) if rest_frame else wave_array
        corr = reference_corr(wave_law, e_bv_test, law)

        assert np.allclose(physical_flux(spec_out), flux_array * corr, rtol=1e-8, equal_nan=True)
        assert np.allclose(spec_out.err_flux.data * norm_flux, err_array * corr, rtol=1e-8, equal_nan=True)

        # Original mask plus the pixels where the law is not defined
        assert np.array_equal(spec_out.flux.mask, pixel_mask | ~np.isfinite(corr))

        return

    def test_default_frame_is_observed(self):

        spec_default = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test)
        spec_obs = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test, rest_frame=False)
        spec_rest = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test, rest_frame=True)

        assert np.allclose(spec_default.flux.data, spec_obs.flux.data, equal_nan=True)
        assert not np.allclose(spec_default.flux.data, spec_rest.flux.data, equal_nan=True)

        # Both frames are the same if there is no redshift
        spec_z0 = lime.Spectrum(wave_array, flux_array, err_array, redshift=0, norm_flux=norm_flux, pixel_mask=pixel_mask)
        out_obs = spectrum_reddening_correction(spec_z0, law='CCM89', E_BV=e_bv_test, rest_frame=False)
        out_rest = spectrum_reddening_correction(spec_z0, law='CCM89', E_BV=e_bv_test, rest_frame=True)
        assert np.allclose(out_obs.flux.data, out_rest.flux.data, equal_nan=True)

        return

    @pytest.mark.parametrize('law', EXT_LAWS)
    def test_equivalent_extinction_inputs(self, law):

        r_v = effective_rv(law)
        klam_Hbeta = reference_klam(HBETA_WAVE, law)

        a_v = e_bv_test * r_v
        c_hbeta = e_bv_test * klam_Hbeta / 2.5

        out_ebv = spectrum_reddening_correction(spec, law=law, E_BV=e_bv_test)
        out_av = spectrum_reddening_correction(spec, law=law, A_V=a_v)
        out_chbeta = spectrum_reddening_correction(spec, law=law, cHbeta=c_hbeta)

        assert np.allclose(out_ebv.flux.data, out_av.flux.data, rtol=1e-8, equal_nan=True)
        assert np.allclose(out_ebv.flux.data, out_chbeta.flux.data, rtol=1e-8, equal_nan=True)
        assert np.array_equal(out_ebv.flux.mask, out_av.flux.mask)
        assert np.array_equal(out_ebv.flux.mask, out_chbeta.flux.mask)

        return

    def test_cHbeta_gives_Hbeta_correction(self):

        # At the Hbeta wavelength the correction must be 10^cHbeta
        c_hbeta = 0.4
        spec_out = spectrum_reddening_correction(spec, law='CCM89', cHbeta=c_hbeta, rest_frame=True)

        idx = np.argmin(np.abs(spec.wave_rest.data - HBETA_WAVE))
        assert not spec_out.flux.mask[idx]
        assert np.isclose(physical_flux(spec_out)[idx] / flux_array[idx], 10 ** c_hbeta, rtol=5e-3)

        return

    def test_pyneb_correction_factor(self):

        # Compare against the correction factor from PyNeb (inside the CCM89 optical-NIR validity range)
        spec_out = spectrum_reddening_correction(spec, law='CCM89', R_V=3.1, E_BV=e_bv_test)
        corr_pn = pn.RedCorr(E_BV=e_bv_test, R_V=3.1, law='CCM89').getCorr(wave_array)

        idcs = (wave_array > 3200) & (wave_array < 9000) & ~pixel_mask
        assert np.allclose(physical_flux(spec_out)[idcs], flux_array[idcs] * corr_pn[idcs], rtol=1e-6)

        return

    def test_default_R_V(self):

        for law in EXT_LAWS:
            out_default = spectrum_reddening_correction(spec, law=law, E_BV=e_bv_test)
            out_rv = spectrum_reddening_correction(spec, law=law, R_V=DEFAULT_RV[law], E_BV=e_bv_test)
            assert np.allclose(out_default.flux.data, out_rv.flux.data, equal_nan=True)

        return

    def test_R_V_changes_the_correction(self):

        out_31 = spectrum_reddening_correction(spec, law='CCM89', R_V=3.1, E_BV=e_bv_test)
        out_50 = spectrum_reddening_correction(spec, law='CCM89', R_V=5.0, E_BV=e_bv_test)

        corr = reference_corr(wave_array, e_bv_test, 'CCM89', r_v=5.0)
        assert np.allclose(physical_flux(out_50), flux_array * corr, rtol=1e-8, equal_nan=True)
        assert not np.allclose(out_31.flux.data, out_50.flux.data, equal_nan=True)

        # A_V uses the requested R_V: A_V = R_V * E(B-V)
        out_av = spectrum_reddening_correction(spec, law='CCM89', R_V=5.0, A_V=5.0 * e_bv_test)
        assert np.allclose(out_av.flux.data, out_50.flux.data, rtol=1e-8, equal_nan=True)

        return

    @pytest.mark.parametrize('law', sorted(FIXED_RV_LAWS))
    def test_fixed_R_V_laws(self, law):

        # The requested R_V is ignored for these laws
        out_default = spectrum_reddening_correction(spec, law=law, E_BV=e_bv_test)
        out_rv = spectrum_reddening_correction(spec, law=law, R_V=5.0, E_BV=e_bv_test)
        assert np.allclose(out_default.flux.data, out_rv.flux.data, equal_nan=True)

        # ... and also in the A_V conversion
        out_av = spectrum_reddening_correction(spec, law=law, R_V=5.0, A_V=e_bv_test * DEFAULT_RV[law])
        assert np.allclose(out_default.flux.data, out_av.flux.data, rtol=1e-8, equal_nan=True)

        return

    def test_zero_extinction(self):

        spec_out = spectrum_reddening_correction(spec, law='CCM89', E_BV=0)

        assert np.allclose(physical_flux(spec_out), flux_array, equal_nan=True)
        assert np.allclose(spec_out.err_flux.data * norm_flux, err_array, equal_nan=True)
        assert np.array_equal(spec_out.flux.mask, pixel_mask)

        return

    def test_correction_increases_flux_towards_the_blue(self):

        spec_out = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test)

        idcs_good, ratio = correction_ratio(spec_out, spec)
        wave_good = spec.wave_rest.data[idcs_good]
        idcs_range = (wave_good > 3500) & (wave_good < 9000)

        assert np.all(ratio[idcs_range] > 1)

        # Monotonic: the bluer the pixel the larger the correction
        order = np.argsort(wave_good[idcs_range])
        assert np.all(np.diff(ratio[idcs_range][order]) < 0)

        return

    def test_uncertainty_scales_like_flux(self):

        spec_out = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test)

        # Flux correction at the good pixels (nan elsewhere)
        idcs_good, ratio_flux = correction_ratio(spec_out, spec)
        ratio_full = np.full(wave_array.size, np.nan)
        ratio_full[idcs_good] = ratio_flux

        # The uncertainty must be scaled by the same factor
        err_in, err_out = spec.err_flux.data * norm_flux, spec_out.err_flux.data * norm_flux
        idcs_err = idcs_good & np.isfinite(err_in) & (err_in != 0)

        assert idcs_err.sum() > 0
        assert np.allclose(err_out[idcs_err] / err_in[idcs_err], ratio_full[idcs_err], rtol=1e-8)

        return

    def test_spectrum_without_uncertainty(self):

        spec_no_err = lime.Spectrum(wave_array, flux_array, redshift=redshift, norm_flux=norm_flux, pixel_mask=pixel_mask)
        spec_out = spectrum_reddening_correction(spec_no_err, law='CCM89', E_BV=e_bv_test)
        spec_ref = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test)

        assert spec_out.err_flux is None
        assert np.allclose(spec_out.flux.data, spec_ref.flux.data, equal_nan=True)

        return

    def test_normalization_preserved(self):

        for norm in (1e-17, 1e-16, 1):
            spec_norm = lime.Spectrum(wave_array, flux_array, err_array, redshift=redshift, norm_flux=norm,
                                      pixel_mask=pixel_mask)
            spec_out = spectrum_reddening_correction(spec_norm, law='CCM89', E_BV=e_bv_test)
            corr = reference_corr(wave_array, e_bv_test, 'CCM89')

            assert spec_out.norm_flux == norm
            assert np.allclose(physical_flux(spec_out), flux_array * corr, rtol=1e-8, equal_nan=True)

        return

    def test_wavelength_units(self):

        # The law must be evaluated in Angstroms whatever the spectrum wavelength units
        spec_nm = lime.Spectrum(wave_array, flux_array, err_array, redshift=redshift, norm_flux=norm_flux,
                                pixel_mask=pixel_mask)
        spec_nm.unit_conversion(wave_units_out='nm')
        assert np.allclose(spec_nm.wave.data, wave_array / 10)

        out_nm = spectrum_reddening_correction(spec_nm, law='CCM89', E_BV=e_bv_test)
        out_AA = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test)

        idcs_nm, ratio_nm = correction_ratio(out_nm, spec_nm)
        idcs_AA, ratio_AA = correction_ratio(out_AA, spec)

        assert np.array_equal(idcs_nm, idcs_AA)
        assert np.allclose(ratio_nm, ratio_AA, rtol=1e-8)
        assert out_nm.units_wave == spec_nm.units_wave

        return

    @pytest.mark.parametrize('kwargs', [{},
                                        {'A_V': 1.0, 'E_BV': 0.3},
                                        {'E_BV': 0.3, 'cHbeta': 0.4},
                                        {'A_V': 1.0, 'cHbeta': 0.4},
                                        {'A_V': 1.0, 'E_BV': 0.3, 'cHbeta': 0.4}])
    def test_one_extinction_parameter_required(self, kwargs):

        with pytest.raises(SpecSyError):
            spectrum_reddening_correction(spec, law='CCM89', **kwargs)

        return

    def test_unknown_law(self):

        with pytest.raises(ValueError):
            spectrum_reddening_correction(spec, law='not_a_law', E_BV=e_bv_test)

        return

    def test_plot_results(self, tmp_path, monkeypatch):

        # Do not open a window and record the arguments sent to the input spectrum plot
        plot_calls = []
        original_plot = type(spec.plot).spectrum

        def record_plot(self, *args, **kwargs):
            plot_calls.append(kwargs)
            return original_plot(self, *args, **kwargs)

        monkeypatch.setattr(type(spec.plot), 'show', lambda self, *args, **kwargs: None)
        monkeypatch.setattr(type(spec.plot), 'spectrum', record_plot)

        fname = tmp_path / 'reddening_comparison.png'
        spec_out = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test, plot_results=True,
                                                 fname=fname, log_scale=False)

        # The output is the same as without the plot
        spec_ref = spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test)
        assert np.allclose(spec_out.flux.data, spec_ref.flux.data, equal_nan=True)

        # "fname" is not forwarded to the plot, the user arguments are, the defaults fill the rest
        assert len(plot_calls) == 1
        assert 'fname' not in plot_calls[0]
        assert plot_calls[0]['log_scale'] is False
        assert plot_calls[0]['rest_frame'] is False

        # The figure was saved and has the corrected spectrum on top of the input one
        assert fname.is_file()
        assert 'Reddening corrected' in [text.get_text() for text in spec.plot.ax.get_legend().get_texts()]

        return

    def test_no_plot_by_default(self, monkeypatch):

        plot_calls = []
        monkeypatch.setattr(type(spec.plot), 'spectrum', lambda self, *args, **kwargs: plot_calls.append(kwargs))

        spectrum_reddening_correction(spec, law='CCM89', E_BV=e_bv_test)
        assert len(plot_calls) == 0

        return