from typing import Any, List, Union

def assert_same_length(
    *items: Union[List[Any], str],
    message : str = "Items should have the same length"
):
    assert len(items) > 0, "Expected at least one list"

    result = len(items[0])

    for item in items:
        assert len(item) == result, message

    return result