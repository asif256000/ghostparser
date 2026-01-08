from ghostparser.parser import generate_response


def test_generate_response():
    """Test the generate_response function."""
    user_input = "Hello, Ghostparser!"
    expected_response = "Response to: Hello, Ghostparser!"
    assert generate_response(user_input) == expected_response
