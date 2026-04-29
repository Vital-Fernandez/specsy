from specsy.models.chemistry import DirectMethod


class ModelManager:

    def __init__(self, lines_df, ion_struct):

        self.direct_method = DirectMethod(lines_df=lines_df, ion_struct=ion_struct)

        return