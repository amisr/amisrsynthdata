Contributing
============

Bug reports, feature suggestions and other contributions are greatly
appreciated! While I can't promise to implement everything, I will always try
to respond in a timely manner.  Please ensure all bug reports, feature
requests, and contribution discussions abide by the
`Code of Conduct <https://amisrsynthdata.readthedocs.io/en/latest/conduct.html#>`_.

Short version
-------------

* Submit bug reports and feature requests at
  `GitHub <https://github.com/amisr/amisrsynthdata/issues>`_
* Make pull requests to the ``develop`` branch

Bug reports
-----------

When `reporting a bug <https://github.com/amisr/amisrsynthdata/issues>`_ please
include:

* Detailed steps to reproduce the bug, including a minimal working example
  where possible
* Full traceback of the error message (when applicable)
* Any details about your local setup that might be helpful in troubleshooting

Feature requests and feedback
-----------------------------

The best way to send feedback is to file an issue at
`GitHub <https://github.com/amisr/amisrsynthdata/issues>`_.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions
  are welcome :)

Development
-----------

To add new features or fix existing bugs in amisrsynthdata, and contribute
them back to the main repository:

1. Set up amisrsynthdata locally for development by following the `development 
   installation <https://amisrsynthdata.readthedocs.io/en/latest/installation.html##development>`_
   instructions. 

3. Create a branch for local development based off of the ``develop`` branch

.. code-block::

    git checkout -b name-of-new-feature origin/develop

4. Make changes locally. Add tests for bugs and new features in the relevant
   test file in the ``tests`` directory. The tests are run with ``pytest``
   and can be written as normal functions (starting with ``test_``)
   containing a standard ``assert`` statement for testing output.

4. Run ``pytest`` locally

.. code-block::

    pytest tests/

5. Commit your changes and push your branch to GitHub

.. code-block::

    git add .
    git commit -m "Brief description of your changes"
    git push origin name-of-new-feature


6. Submit a pull request through the GitHub website. Pull requests should be
   made to the ``develop`` branch. The continuous integration (CI) testing
   servers will automatically test the whole codebase, including your changes,
   for multiple versions of Python.

Pull Request Guidelines
^^^^^^^^^^^^^^^^^^^^^^^

If you need some code review or feedback while you're developing the code, just
make a pull request.

For merging, you should:

1. Include passing tests for your changes
2. Update/add documentation if relevant
3. Include a detailed description of the changes and their impacts in the pull
   request description

Style Guidelines
^^^^^^^^^^^^^^^^

In general, `amisrsynthdata` follows PEP8 and numpydoc guidelines.  PyTest is
used to run the unit and integration tests, flake8 checks for style, and
sphinx-build performs documentation tests.  However, there are certain
additional style elements that have been settled on to ensure the project
maintains a consistent coding style:

- Line breaks should occur before a binary operator (ignoring flake8 W503)
- Preferably break long lines on open parentheses instead of using backslashes
- Use no more than 80 characters per line
- Several dependent packages have common nicknames, including:

  * ``import datetime as dt``
  * ``import numpy as np``

- Provide tests with informative failure statements and descriptive, one-line
  docstrings.

Add Ionospheric State Function
------------------------------

The package was deliberately designed to make it straight forwards for users to expand the libary of ionospheric state functions by adding new ones.  Instructions for how to do this are available in the `Create New State Functions <https://amisrsynthdata.readthedocs.io/en/latest/ionostate.html#create-new-state-functions>`_ section.
