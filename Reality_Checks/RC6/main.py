from suggested_activities import *


def main():
    # Display available options to the user
    choice = input("Indicate which Suggested Activity you would like to run (1-6, 0 to run all):")

    # Execute a function dependent on user input
    if choice == "1":
        sa1()
    elif choice == "2":
        print("Suggested Activity 2")
        sa2()
    elif choice == "3":
        print("Suggested Activity 3")
        sa3()
    elif choice == "4":
        print("Suggested Activity 4")
        sa4()
    elif choice == "5":
        print("Suggested Activity 5")
        sa5()
    elif choice == "6":
        print("Suggested Activity 6")
        sa6()
    elif choice == "0":
        print("All Activities")
        sa1()
        sa2()
        sa3()
        sa4()
        sa5()
        sa6()
    else:
        print("Invalid input.")

    return 0


# While not required, it is considered good practice to have
# a main function and use this syntax to call it.
if __name__ == "__main__":
    main()
