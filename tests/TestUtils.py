import os


def get_tests_directory():
    """discover and return the absolute path of the root directory"""
    testdir = os.getcwd()
    if not testdir.endswith("tests"):
        testdir = os.path.join(testdir, "tests")
        if not os.path.exists(testdir):
            raise ValueError("can not find test directory")

    return testdir
