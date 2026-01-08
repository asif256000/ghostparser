"""Entry point for the ghostparser module."""

from ghostparser.parser import generate_response


def main():
    """Main entry point for the ghostparser module."""
    user_input = input("Enter your input: ")
    response = generate_response(user_input)
    print(response)


if __name__ == "__main__":
    main()
