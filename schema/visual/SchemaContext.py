class SchemaContext():

    __instances = 0

    def __init__(self) -> None:
        SchemaContext.__instances += 1

    def raise_error_message(self, error):
        print(error)

    @staticmethod
    def is_loaded() -> bool:
        return SchemaContext.__instances > 0

    @property
    def clustal(self):
        return "clustalo"
