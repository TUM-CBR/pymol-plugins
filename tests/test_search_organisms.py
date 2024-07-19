from tests.common import *

class TestSearchOrganisms:
    def test_search_organisms(self, pymol_fixture):
        from time import sleep
        import threading
        print("Test thread: ", threading.current_thread().ident)
        sleep(2)
        import pymol
        pymol.cmd.fetch('1aon')
        pymol.cmd.color('blue', '1aon')
        sleep(2)

