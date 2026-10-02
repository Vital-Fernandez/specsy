import pytest
import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import matplotlib
import logging

# matplotlib.use('Agg')

from specsy.models.ssp import Binary, StellarBinaries, to_dataframe
from specsy.io import SpecSyError, specsy_cfg
from matplotlib import pyplot as plt

SOURCE, LIBRARY, IMF = 'bpass', 'AP', 'imf135_300'
AGES = np.round(np.arange(6.0, 11.01, 0.1), 1)
METALS = np.array([1e-5, 1e-4, 1e-3, 0.002, 0.003, 0.004, 0.006, 0.008, 0.010, 0.014, 0.020, 0.030, 0.040])
WAVE = np.arange(912.0, 3001.0, 1.0)

script_dir = Path(__file__).resolve().parent
bpass_pname = script_dir / 'baseline' / 'SESAMME_BAPSS-AP_v3.fits'
bpass_grid = StellarBinaries.from_fits(bpass_pname) if bpass_pname.exists() else None

# print(bpass_grid)
# bpass_grid.plot_age_metallicity()
# wave, flux = bpass_grid.get_spectrum(8.05, 0.012, interpolate=True, plot_results=False)


def synth_flux(age, met, offset=0.0):
    return age * 1e3 + met * 1e5 + 1e-3 * WAVE + offset


def make_binaries(loaded=True, library=LIBRARY, offset=0.0, imf=IMF):
    return [Binary(source=SOURCE,
                   library=library,
                   metallicity=float(z),
                   alpha=0.0, age=float(a),
                   fpath=f'synthetic_{z}_{a}.dat',
                   imf=imf,
                   flux=synth_flux(a, z, offset) if loaded else None,
                   wavelength=WAVE.copy() if loaded else None) for z in METALS for a in AGES]


@pytest.fixture
def sb():
    return StellarBinaries(make_binaries())


def _assert_reference_grid(sb, n_loaded=663):
    assert len(sb.frame) == 663
    np.testing.assert_allclose(sb.ages, AGES)
    np.testing.assert_allclose(sb.metallicities, METALS)
    assert sb.source == SOURCE and sb.library == LIBRARY
    assert sb.uniform_wmin == pytest.approx(912.0)
    assert sb.uniform_wmax == pytest.approx(3000.0)
    assert sb.uniform_deltalambda == pytest.approx(1.0)
    assert sb.uniform_dispersion and sb.is_loaded == (n_loaded == 663)
    assert list(sb.wave_series.index) == [0] and sb.wave_series[0].size == WAVE.size
    assert f'Total stellar spectra : 663 ({n_loaded} loaded)' in repr(sb)
    return


# --------------- Binary Class tests
def test_binary_load_state():
    b = make_binaries()[0]
    assert b.is_loaded
    b.clear_spectrum()
    assert not b.is_loaded and b.wavelength is None
    return


# --------------- Stellar binaries declaration tests

def test_reference_grid(sb):
    _assert_reference_grid(sb)
    return


def test_frame_structure(sb):
    assert list(sb.frame.columns) == specsy_cfg['stellar']['ssp_params']
    assert sb.frame.index.equals(pd.RangeIndex(663))
    assert sb.flux_series.index.equals(sb.frame.index)
    np.testing.assert_allclose(sb.flux_series[0], synth_flux(AGES[0], METALS[0]))
    return


def test_to_dataframe(sb):
    pd.testing.assert_frame_equal(to_dataframe(make_binaries()), sb.frame)
    df = to_dataframe(make_binaries(), columns=['age', 'metallicity'])
    assert list(df.columns) == ['age', 'metallicity'] and len(df) == 663
    return


def test_save_frame(sb, tmp_path):
    sb.save_frame(tmp_path / 'frame.txt')
    assert (tmp_path / 'frame.txt').exists()
    return


def test_unloaded_cube():
    sb = StellarBinaries(make_binaries(loaded=False))
    assert len(sb.frame) == 663 and sb.ages.size == 51 and sb.metallicities.size == 13
    assert sb.flux_series is None and sb.wave_series is None and not sb.is_loaded
    assert sb.uniform_wmin is None and not sb.uniform_dispersion
    assert sb.source == SOURCE and sb.library == LIBRARY
    assert 'Total stellar spectra : 663 (0 loaded)' in repr(sb)
    return


def test_partially_loaded_cube(caplog):
    binaries = make_binaries()
    binaries[0].clear_spectrum()
    with caplog.at_level(logging.WARNING):
        sb = StellarBinaries(binaries)
    assert not sb.is_loaded and '1 of the 663 stellar spectra are not loaded' in caplog.text
    assert sb.get_binaries(age=6.0, metallicity=1e-5)[0].flux is None
    np.testing.assert_allclose(sb.get_spectrum(6.1, 1e-5)[1], synth_flux(6.1, 1e-5))
    return


def test_duplicated_rows_warns(caplog):
    with caplog.at_level(logging.WARNING):
        StellarBinaries(make_binaries() + make_binaries()[:3])
    assert '3 of the 666 rows' in caplog.text
    return


def test_nonuniform_grid_warns(caplog):
    binaries = make_binaries()
    binaries[0].wavelength = binaries[0].wavelength + 0.5
    with caplog.at_level(logging.WARNING):
        sb = StellarBinaries(binaries)
    assert sb.uniform_wmin is None and sb.uniform_wmax is None
    assert sb.uniform_deltalambda == pytest.approx(1.0)
    assert not sb.uniform_dispersion and sb.wave_series.index.equals(sb.frame.index)
    assert 'not uniform' in caplog.text and 'wmin' in caplog.text
    return


def test_multiple_libraries_warns(caplog):
    with caplog.at_level(logging.WARNING):
        sb = StellarBinaries(make_binaries() + make_binaries(library='C3K'))
    assert sb.library == ['AP', 'C3K'] and sb.source == SOURCE
    assert 'multiple libraries' in caplog.text
    return


# Selection and indexing
def test_get_binaries(sb):
    assert len(sb.get_binaries()) == 663
    assert len(sb.get_binaries(source='bpass')) == 663
    assert len(sb.get_binaries(source='bpass', library='AP', imf=IMF)) == 663
    assert len(sb.get_binaries(age=8.0)) == 13
    assert len(sb.get_binaries(age=[7.0, 8.0])) == 26

    b = sb.get_binaries(age=8.0, metallicity=0.014)[0]
    assert isinstance(b, Binary) and b.library == LIBRARY and b.imf == IMF and b.low_imf_exp is None
    np.testing.assert_array_equal(b.wavelength, WAVE)
    np.testing.assert_allclose(b.flux, synth_flux(8.0, 0.014))
    return


@pytest.mark.parametrize('kwargs, n', [({'age': 6.1}, 13),
                                       ({'metallicity': 1e-5}, 51),
                                       ({'metallicity': 1e-4}, 51),
                                       ({'age': [6.0, 6.1], 'metallicity': 0.014}, 2),
                                       ({'library': 'AP', 'age': 11.0}, 13)])
def test_exact_selection(sb, kwargs, n):
    assert len(sb.get_binaries(**kwargs)) == n
    return


@pytest.mark.parametrize('kwargs, match', [({'age': 6.1 + 1e-7}, 'not found'),  # exact match, no tolerance
                                           ({'age': 12.0}, 'not found'),
                                           ({'metallicity': 1.1e-5}, 'not found'),
                                           ({'source': 'BPASS'}, 'not found'),  # case-sensitive
                                           ({'colour': 'red'}, 'Unrecognised')])
def test_selection_errors(sb, kwargs, match):
    with pytest.raises(SpecSyError, match=match):
        sb.get_binaries(**kwargs)
    return


def test_selection_missing_combination():
    # C3K only covers the youngest age: each value exists, but not together
    binaries = make_binaries() + [b for b in make_binaries(library='C3K') if b.age == 6.0]
    sb = StellarBinaries(binaries)
    with pytest.raises(SpecSyError, match='No stellar spectra match'):
        sb.get_binaries(library='C3K', age=7.0)
    return


def test_get_spectrum_matches_fits():

    age, met = 8.0, 0.014

    with fits.open(bpass_pname) as hdul:
        hdu = next(h for h in hdul[1:] if np.isclose(h.header['METAL'], met))
        wave_ref = np.array(hdu.data['WL'])
        flux_ref = np.array(hdu.data[str(age)])

    wave, flux = bpass_grid.get_spectrum(age, met)

    np.testing.assert_array_equal(wave, wave_ref)
    np.testing.assert_allclose(flux, flux_ref)

    return


def test_get_spectrum_errors():
    with pytest.raises(SpecSyError, match='not found'):
        StellarBinaries(make_binaries()).get_spectrum(12.0, 0.014)

    sb2 = StellarBinaries(make_binaries() + make_binaries(library='C3K', offset=1.0))
    with pytest.raises(SpecSyError, match='Multiple stellar spectra'):
        sb2.get_spectrum(8.0, 0.014)
    np.testing.assert_allclose(sb2.get_spectrum(8.0, 0.014, library='C3K')[1], synth_flux(8.0, 0.014, 1.0))

    with pytest.raises(SpecSyError, match='not loaded'):
        StellarBinaries(make_binaries(loaded=False)).get_spectrum(8.0, 0.014)
    return


@pytest.mark.parametrize('limits, n_mets', [(None, 13), ((None, None), 13),
                                            ((0.001, 0.02), 7),  # bounds are exclusive
                                            ((None, 0.02), 10), ((0.001, None), 10)])
def test_limits_mask(sb, limits, n_mets):
    assert sb._limits_mask('metallicity', limits).sum() == n_mets * 51
    return


# --------------- FITS export / import tests
@pytest.mark.parametrize('kwargs, match', [({'met_limits': (0.02,)}, 'tuple'),
                                           ({'met_limits': (0.02, 0.001)}, 'lower value'),
                                           ({'age_limits': (8.0, 8.0)}, 'lower value'),
                                           ({'age_limits': (11.0, None)}, 'No spectra left'),
                                           ({'libraries': 'C3K'}, 'No valid libraries'),
                                           ({'pixel_width': 2.0, 'pixel_number': 100}, 'mutually exclusive'),
                                           ({'idcs_ssp': [0, 1000]}, 'not in the frame')])
def test_to_fits_errors(sb, tmp_path, kwargs, match):
    with pytest.raises(SpecSyError, match=match):
        sb.to_fits(tmp_path / 'cube.fits', **kwargs)
    return


def test_to_fits_missing_folder_and_unloaded(sb, tmp_path):
    with pytest.raises(SpecSyError, match='Output folder not found'):
        sb.to_fits(tmp_path / 'missing' / 'cube.fits')
    with pytest.raises(SpecSyError, match='not loaded'):
        StellarBinaries(make_binaries(loaded=False)).to_fits(tmp_path / 'cube.fits')
    return


def test_fits_save_load(sb, tmp_path):
    fname = tmp_path / 'bpass_AP.fits'
    sb.to_fits(fname)

    with fits.open(fname) as hdul:
        assert len(hdul) == 1 + 13
        hdr = hdul[0].header
        assert hdr['SOURCE'] == 'bpass' and hdr['LIBS'] == 'AP' and hdr['IMF'] == IMF
        assert hdr['MET_MIN'] == pytest.approx(1e-5) and hdr['MET_MAX'] == pytest.approx(0.04)
        assert hdr['AGE_MIN'] == pytest.approx(6.0) and hdr['AGE_MAX'] == pytest.approx(11.0)
        assert len(hdul[1].columns) == 1 + 51
        assert [h.header['METAL'] for h in hdul[1:]] == pytest.approx(list(METALS))
        assert all(h.header['LIBRARY'] == LIBRARY for h in hdul[1:])
        assert not any(h.name.endswith(f'_{LIBRARY}') for h in hdul[1:])  # single library: no suffix

    sb2 = StellarBinaries.from_fits(fname)
    _assert_reference_grid(sb2)
    for age, met in [(6.0, 1e-5), (8.0, 0.014), (11.0, 0.04)]:
        np.testing.assert_allclose(sb2.get_spectrum(age, met, imf=IMF)[1], sb.get_spectrum(age, met)[1])

    return


def test_to_fits_limits_and_crop(sb, tmp_path):
    fname = tmp_path / 'subset.fits'
    sb.to_fits(fname, met_limits=(0.001, 0.02), age_limits=(7.0, 9.0), wmin=1200, wmax=2000)

    with fits.open(fname) as hdr:
        assert hdr[0].header['MET_MIN'] == pytest.approx(0.002)
        assert hdr[0].header['AGE_MIN'] == pytest.approx(7.1)

    sb2 = StellarBinaries.from_fits(fname)
    np.testing.assert_allclose(sb2.metallicities, [0.002, 0.003, 0.004, 0.006, 0.008, 0.010, 0.014])
    np.testing.assert_allclose(sb2.ages, AGES[(AGES > 7.0) & (AGES < 9.0)])
    assert sb2.uniform_wmin == 1200.0 and sb2.uniform_wmax == 2000.0
    return


def test_fits_load_save_two_libraries(tmp_path):
    sb = StellarBinaries(make_binaries() + make_binaries(library='C3K', offset=1.0))
    fname = tmp_path / 'two_libs.fits'
    sb.to_fits(fname)

    with fits.open(fname) as hdul:
        assert len(hdul) == 1 + 2 * 13
        assert {h.header['LIBRARY'] for h in hdul[1:]} == {'AP', 'C3K'}
        assert all(h.name.endswith(f'_{h.header["LIBRARY"]}') for h in hdul[1:])

    sb2 = StellarBinaries.from_fits(fname)
    assert len(sb2.frame) == 2 * 663 and sb2.library == ['AP', 'C3K']
    np.testing.assert_allclose(sb2.get_spectrum(8.0, 0.014, library='AP')[1], synth_flux(8.0, 0.014))
    np.testing.assert_allclose(sb2.get_spectrum(8.0, 0.014, library='C3K')[1], synth_flux(8.0, 0.014, 1.0))
    return


def test_to_fits_same_age_raises(tmp_path):
    # Two IMFs within one library would overwrite each other's age columns
    sb = StellarBinaries(make_binaries() + make_binaries(imf='imf100_300', offset=2.0))
    with pytest.raises(SpecSyError, match='share the same age'):
        sb.to_fits(tmp_path / 'two_imfs.fits')
    return


def test_to_fits_resampled(sb, tmp_path):
    fname = tmp_path / 'resampled.fits'
    sb.to_fits(fname, wmin=1200, wmax=2000, pixel_width=2.0)
    sb2 = StellarBinaries.from_fits(fname)
    assert sb2.uniform_deltalambda == pytest.approx(2.0) and len(sb2.frame) == 663
    return


def test_to_fits_idcs_ssp(sb, tmp_path):
    fname = tmp_path / 'z014.fits'
    sb.to_fits(fname, idcs_ssp=sb.frame.index[sb.frame['metallicity'] == 0.014])
    with fits.open(fname) as hdul:
        assert len(hdul) == 1 + 1 and len(hdul[1].columns) == 1 + 51
    return


@pytest.mark.skipif(bpass_grid is None, reason=f'Reference BPASS cube not found: {bpass_pname}')
def test_real_bpass_cube():
    _assert_reference_grid(bpass_grid)
    return


@pytest.mark.skipif(bpass_grid is None, reason=f'Reference BPASS cube not found: {bpass_pname}')
def test_to_fits_to_stellar_binaries_roundtrip(tmp_path):

    """
    to_fits and to_stellar_binaries take the same selection: their outputs should hold the same spectra,
    and both should still match the original reference cube.
    """

    idcs_ssp = bpass_grid.frame.index[bpass_grid.frame['metallicity'].isin([1e-5, 0.008, 0.02])]

    fname = tmp_path / 'subset.fits'
    bpass_grid.to_fits(fname, idcs_ssp=idcs_ssp)
    from_file = StellarBinaries.from_fits(fname)

    from_memory = bpass_grid.to_stellar_binaries(idcs_ssp=idcs_ssp)

    assert len(from_file.frame) == len(from_memory.frame) == len(idcs_ssp)
    np.testing.assert_allclose(from_file.metallicities, from_memory.metallicities)
    np.testing.assert_allclose(from_file.ages, from_memory.ages)

    for age, met in [(bpass_grid.ages[0], 1e-5), (bpass_grid.ages[-1], 0.008), (bpass_grid.ages[len(bpass_grid.ages) // 2], 0.02)]:
        wave_ref, flux_ref = bpass_grid.get_spectrum(age, met)
        wave_file, flux_file = from_file.get_spectrum(age, met)
        wave_mem, flux_mem = from_memory.get_spectrum(age, met)

        np.testing.assert_allclose(wave_file, wave_ref)
        np.testing.assert_allclose(flux_file, flux_ref)
        np.testing.assert_allclose(wave_mem, wave_ref)
        np.testing.assert_allclose(flux_mem, flux_ref)

    return


# --------------- Spectra operation tests

def test_to_stellar_binaries_copy(sb):
    sb2 = sb.to_stellar_binaries()
    assert sb2 is not sb
    _assert_reference_grid(sb2)
    np.testing.assert_allclose(sb2.get_spectrum(8.0, 0.014)[1], synth_flux(8.0, 0.014))
    return


def test_to_stellar_binaries_resampling(sb):
    sb2 = sb.to_stellar_binaries(wmin=1200, wmax=2000, pixel_width=2.0)

    assert sb2 is not sb and len(sb2.frame) == 663
    assert sb2.source == SOURCE and sb2.library == LIBRARY
    assert sb2.uniform_deltalambda == pytest.approx(2.0)
    assert sb2.uniform_dispersion and list(sb2.wave_series.index) == [0]
    assert sb.uniform_deltalambda == pytest.approx(1.0)  # original untouched

    wave, flux = sb2.get_spectrum(8.0, 0.014)
    assert np.all(np.isfinite(flux))
    assert wave.min() >= 1199.0 and wave.max() <= 2001.0
    ref = synth_flux(8.0, 0.014)[(WAVE >= 1200) & (WAVE <= 2000)]
    assert np.mean(flux) == pytest.approx(np.mean(ref), rel=1e-3)  # flux density level preserved

    return


def test_to_stellar_binaries_selection(sb):
    sb2 = sb.to_stellar_binaries(met_limits=(0.001, 0.02), age_limits=(7.0, 9.0), wmin=1200, wmax=2000)
    assert len(sb2.frame) == 7 * 19 and sb2.frame.index.equals(pd.RangeIndex(7 * 19))
    assert sb2.uniform_wmin == 1200.0 and sb2.uniform_wmax == 2000.0 and sb2.uniform_deltalambda == 1.0
    return


def test_to_stellar_binaries_library_resets_frame():
    sb = StellarBinaries(make_binaries() + make_binaries(library='C3K', offset=1.0))
    sb2 = sb.to_stellar_binaries(libraries='C3K', pixel_width=2.0)
    assert sb2.library == 'C3K' and sb2.frame.index.equals(pd.RangeIndex(663))
    assert sb2.flux_series.index.equals(sb2.frame.index)
    return


def test_idcs_ssp(sb):
    # Integer frame indices keep the requested order
    sb2 = sb.to_stellar_binaries(idcs_ssp=[52, 0, 1])
    assert list(sb2.frame['age']) == [AGES[1], AGES[0], AGES[1]] and list(sb2.frame.index) == [0, 1, 2]

    # Boolean mask from the frame
    sb3 = sb.to_stellar_binaries(idcs_ssp=sb.frame['age'] > 10.0)
    assert len(sb3.frame) == 10 * 13
    return


def test_idcs_ssp_precedence(sb, caplog):
    with caplog.at_level(logging.WARNING):
        sb2 = sb.to_stellar_binaries(idcs_ssp=[0, 1], libraries='AP', age_limits=(9.0, 10.0))
    assert len(sb2.frame) == 2 and 'takes precedence' in caplog.text
    return


@pytest.mark.parametrize('idcs_ssp, match', [([0, 663], 'not in the frame'),
                                             (np.ones(10, dtype=bool), 'one entry per frame row'),
                                             ([0.5, 1.0], 'integer frame indices'),
                                             ([], 'No spectra left')])
def test_idcs_ssp_errors(sb, idcs_ssp, match):
    with pytest.raises(SpecSyError, match=match):
        sb.to_stellar_binaries(idcs_ssp=idcs_ssp)
    return


def test_to_stellar_binaries_errors(sb):
    with pytest.raises(SpecSyError, match='mutually exclusive'):
        sb.to_stellar_binaries(pixel_width=2.0, pixel_number=500)
    with pytest.raises(SpecSyError, match='Fewer than 2 pixels'):
        sb.to_stellar_binaries(wmin=1500.0, wmax=1500.5)
    with pytest.raises(SpecSyError, match='not loaded'):
        StellarBinaries(make_binaries(loaded=False)).to_stellar_binaries(pixel_width=2.0)
    with pytest.raises(SpecSyError, match='No valid libraries'):
        sb.to_stellar_binaries(libraries='C3K')
    return

# --------------- Source construction tests

def test_from_source_errors(tmp_path):
    with pytest.raises(SpecSyError, match='not recognized'):
        StellarBinaries.from_source('NOT_A_SOURCE', tmp_path)
    valid_source = specsy_cfg['stellar']['source_list'][0]
    with pytest.raises(SpecSyError, match='Folder not found'):
        StellarBinaries.from_source(valid_source, tmp_path / 'missing')
    return


def test_source_list_names():
    assert 'pystarburst99' in specsy_cfg['stellar']['source_list']
    return


# --------------- Plotting tests

def test_plots_smoke(sb, tmp_path):
    sb.plot_age_metallicity(fname=tmp_path / 'grid.png')
    sb.plot_spectra(fname=tmp_path / 'spec.png', metallicity=0.014, age_range=(7.0, 8.0))
    assert (tmp_path / 'grid.png').exists() and (tmp_path / 'spec.png').exists()
    with pytest.raises(SpecSyError, match='No loaded binaries'):
        sb.plot_spectra(fname=tmp_path / 'none.png', age=12.0)
    return


# --------------- _from_frame


def test_from_frame_reload_saved(sb, tmp_path):
    fname = tmp_path / 'frame.csv'
    sb.save_frame(fname)

    reloaded = StellarBinaries.from_frame(fname)
    assert len(reloaded.frame) == 663 and not reloaded.is_loaded and reloaded.flux_series is None
    assert reloaded.source == SOURCE and reloaded.library == LIBRARY
    np.testing.assert_allclose(reloaded.ages, sb.ages)
    np.testing.assert_allclose(reloaded.metallicities, sb.metallicities)
    # NaN/None on the empty IMF-limit columns may not round-trip identically through a text format
    text_cols = ['source', 'library', 'imf', 'fpath']
    pd.testing.assert_frame_equal(reloaded.frame[text_cols], sb.frame[text_cols])
    numeric_cols = ['age', 'metallicity', 'alpha']
    pd.testing.assert_frame_equal(reloaded.frame[numeric_cols], sb.frame[numeric_cols], check_dtype=False)
    return


def test_from_frame_missing_columns():
    frame = to_dataframe(make_binaries()).drop(columns=['age'])
    with pytest.raises(SpecSyError, match='missing the columns'):
        StellarBinaries.from_frame(frame, [None] * len(frame), [None] * len(frame))
    return


def test_from_frame_unknown_columns_warns(caplog):
    frame = to_dataframe(make_binaries())
    frame['extra'] = 1
    with caplog.at_level(logging.WARNING):
        sb2 = StellarBinaries.from_frame(frame, [None] * len(frame), [None] * len(frame))
    assert 'not recognised' in caplog.text and 'extra' not in sb2.frame.columns
    return


def test_from_frame_nonnumeric_index_warns(caplog):
    frame = to_dataframe(make_binaries())
    frame.index = [f'row{i}' for i in range(len(frame))]
    with caplog.at_level(logging.WARNING):
        sb2 = StellarBinaries.from_frame(frame, [None] * len(frame), [None] * len(frame))
    assert 'index is not numeric' in caplog.text
    assert sb2.frame.index.equals(pd.RangeIndex(len(frame)))
    return


def test_from_frame_quiet_suppresses_warnings(caplog):
    frame = to_dataframe(make_binaries())
    frame.index = [f'row{i}' for i in range(len(frame))]
    with caplog.at_level(logging.WARNING):
        StellarBinaries.from_frame(frame, [None] * len(frame), [None] * len(frame), quiet=True)
    assert caplog.text == ''
    return


def test_interpolate_on_node(sb):
    np.testing.assert_array_equal(sb.get_spectrum(8.0, 0.014, interpolate=True)[1], sb.get_spectrum(8.0, 0.014)[1])
    return


def test_interpolate_age(sb):
    # synth_flux is linear in age, so the interpolation is exact
    wave, flux = sb.get_spectrum(8.05, 0.014, interpolate=True)
    np.testing.assert_array_equal(wave, WAVE)
    np.testing.assert_allclose(flux, synth_flux(8.05, 0.014))
    return


def test_interpolate_metallicity_log(sb):
    # Midpoint in log10(Z) between two nodes: equal weights on both
    met = 10 ** ((np.log10(0.010) + np.log10(0.014)) / 2)
    flux = sb.get_spectrum(8.0, met, interpolate=True)[1]
    np.testing.assert_allclose(flux, 0.5 * synth_flux(8.0, 0.010) + 0.5 * synth_flux(8.0, 0.014))
    return


def test_interpolate_with_params():
    sb2 = StellarBinaries(make_binaries() + make_binaries(library='C3K', offset=1.0))
    flux = sb2.get_spectrum(8.05, 0.014, interpolate=True, library='C3K')[1]
    np.testing.assert_allclose(flux, synth_flux(8.05, 0.014, 1.0))
    with pytest.raises(SpecSyError, match='Multiple stellar spectra'):
        sb2.get_spectrum(8.05, 0.014, interpolate=True)
    return


@pytest.mark.parametrize('age, met', [(5.9, 0.014), (11.1, 0.014), (8.0, 1e-6), (8.0, 0.05)])
def test_interpolate_outside_grid(sb, age, met):
    with pytest.raises(SpecSyError, match='outside the grid coverage'):
        sb.get_spectrum(age, met, interpolate=True)
    return


def test_interpolate_missing_node():
    binaries = [b for b in make_binaries() if not (b.age == 8.1 and b.metallicity == 0.014)]
    sb = StellarBinaries(binaries)
    with pytest.raises(SpecSyError, match='is missing'):
        sb.get_spectrum(8.05, 0.014, interpolate=True)
    return


def test_no_interpolation_by_default(sb):
    with pytest.raises(SpecSyError, match='not found'):
        sb.get_spectrum(8.05, 0.014)
    return


def test_interpolate_plot(sb, monkeypatch):
    monkeypatch.setattr(plt, 'show', lambda: None)  # keep the figure open to inspect it
    plt.close('all')

    wave, flux = sb.get_spectrum(8.05, 0.012, interpolate=True, plot_results=True)
    assert len(plt.get_fignums()) == 1

    ax = plt.gcf().axes[0]
    lines = ax.get_lines()
    assert len(lines) == 5  # four grid nodes + the interpolated spectrum

    # The last line is the returned spectrum, the node lines are grid spectra
    np.testing.assert_allclose(lines[-1].get_xdata(), wave)
    np.testing.assert_allclose(lines[-1].get_ydata(), flux)
    np.testing.assert_allclose(lines[0].get_ydata(), sb.get_spectrum(8.0, 0.010)[1])

    # Legend: node weights add up to one and the interpolated entry comes last
    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert labels[-1].startswith('Interpolated')
    node_weights = [float(label.split('w=')[1].rstrip(')')) for label in labels[:-1]]
    assert sum(node_weights) == pytest.approx(1.0, abs=2e-3)

    plt.close('all')
    return