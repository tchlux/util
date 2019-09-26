


if __name__ == "__main__":
    # Modify the path to ensure that the files in this directory are used.
    import os, sys
    cwd = os.path.dirname(os.path.abspath(__file__))
    sys.path = [cwd] + sys.path
    # Tests each module.
    from test_random import test; test()
    from test_stats import test; test()
    from test_system import test; test()
    from test_approximate import test; test()

