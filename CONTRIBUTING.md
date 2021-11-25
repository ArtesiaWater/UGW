# Contributing to UGW: Bouwsteen Oppervlaktewater

We welcome you to contribute code and documentation to UGW! This section
describes how you can contribute to UGW and how the process works of
finally implementing your code. GitHub, where UGW is hosted, also has 
[good tutorials](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests)
to learn how to commit code changes to GitHub open source projects. Let's start!

## 1. Create a GitHub Issue

Before you start you can start a GitHub Issue describing the changes you
propose to make and why these are necessary. This is an easy way to inform
the UGW community in an early stage of any issues that need to be solved
and allows others to help you work out a solution.

## 2. Fork UGW

To start making changes to the original code, you need to make a local copy
of the UGW, called "Forking" in git-language. You can read how to fork a
GitHub repository [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo).

_**Note**:
    Make sure to make changes in the make changes in your local Development
    branch (dev) or start an entirely new branch that branches of the
    dev-branch._

## 3. Write Code

After you forked UGW, you can start making changes to the code or add new
features to it. To ensure high quality code that is easy to read and maintain
we follow [Python PEP8](https://www.python.org/dev/peps/pep-0008/)
standard.

## 4. Document Code

When submitting a new function, method or class, docstrings are required
before the new code will be pulled into the dev branch. Docstrings within
the method or class need to be written in 
[NumPy docformat](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard). For code not in a function method or class it is recommended to add comments to ensure
the code is legible and can be understood by others.

## 6. Create a pull request

Once you have written, tested, and documented your code you can start a pull
request on the development branch (dev) of UGW. Pull requests can only
be submitted to the dev-branch and need to be reviewed by one of the core
developers. After you have create a Pull Request the Core Development Team will
review your code and discuss potential improvements on GitHub before merging
your code into the development branch. After a successful Pull Request your
code will be included in the next release!
