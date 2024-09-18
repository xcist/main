import os


def functional_test_disabled():
    if os.getenv("ENABLE_FUNCTIONAL_TEST", "TRUE") == "FALSE":
        print("Functional Test Disabled")
        return True
    return False
