How to doc
==========

.. toctree::
    :hidden:

    self
    Plugins <dev_plugin>

Just create a doc string for each function, method, class and module or contribute by writing or enhancing these doc strings for existing entities.

If you write a new module, please include it in the API documentation: Add a corresponding RST file in the ``docs/src`` folder and reference the new module in the RST file of the parent module.

Check your documentation by building it locally, look out for new reference errors.

.. rubric:: Build the documentation locally

The documentation is hosted at ReadTheDocs and built automatically on every push to GitHub.

To build the documentation locally:

1. Install the Python packages ``sphinx``, ``furo`` and ``sphinx_autorun``
2. Clone the sugar repository
3. Install sugar, i.e. in editable mode ``cd sugar; pip install -e .``
4. Change to the ``docs`` directory
5. Build the documentation ``sphinx-build -anE . _build``
6. The documentation is now located in the ``docs/_build`` directory
