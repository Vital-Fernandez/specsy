from models.chemistry import DirectMethod


class ReviewInputs:

    def __init__(self, spectrum=None, lines_frame=None, emiss_db=None, R_V=None, stellar_db=None, extinction_law=None):

        return


class Treatments:

    def __init__(self):


        return

    def direct_method(self, line_list):


        return


class Observation:

    def __init__(self, spectrum=None, lines_frame=None, emiss_db=None, R_V=None, stellar_db=None, extinction_law=None,
                 id_name=None):

        # Declare attributes
        self.id_name = id_name

        # Review the data before the measuremetns
        self.data = ReviewInputs(lines_frame, emiss_db, R_V, extinction_law)

        # Treatment objects
        self.fit = Treatments()

        return