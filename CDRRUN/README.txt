README
------

Acknowledging CMIP
==================

The data provided here are derived from the outputs of the `Coupled Model Intercomparison Project <https://www.wcrp-climate.org/wgcm-cmip>`_ data.
This means that you **must** abide by the terms of use of the data, in particular the required acknowledgement statements (see the `CMIP5 terms of use <https://pcmdi.llnl.gov/mips/cmip5/terms-of-use.html>`_ and `CMIP6 terms of use <https://pcmdi.llnl.gov/CMIP6/TermsOfUse/TermsOfUse6-1.html>`_).

To make it easier to do this, we have developed some basic tools which simplify the process of checking model license terms and creating the tables required in publications to cite CMIP data (see the `our tools here <https://netcdf-scm.readthedocs.io/en/latest/usage/using-cmip-data.html>`_).
However, we provide no guarantees that these tools are up to date so all users should double check that they do in fact produce output consistent with the terms of use referenced above (and if there are issues, please raise an issue at netCDF-SCM's `issue tracker <https://gitlab.com/netcdf-scm/netcdf-scm/issues>`_).

In addition, we make no guarantees about the correctness of this data.
We have described the process used to generate this dataset in [paper link will be placed here when available].
The source code required for all calculations is openly available at `<gitlab.com/netcdf-scm/netcdf-scm>`_ (netCDF-SCM releases are archived at `<https://doi.org/10.5281/zenodo.3903372>`_) and `<https://gitlab.com/netcdf-scm/calibration-data>`_ ([zenodo release to come]).

Citation
========

If you use this data, please cite the relevant publication using your format of choice (further formats are available at [paper link will be placed here when available]).

Plain text
++++++++++

Nicholls, Z, Lewis, J, Makin, M, Nattala, U, Zhang, GZ, Mutch, SJ, Tescari, E, Meinshausen, M. Regionally aggregated, stitched and de-drifted CMIP-climate data, processed with netCDF-SCM v2.0.0. Geosci Data J. 2020; submitted.

Bibtex
++++++

@article{nicholls_gdj2020_cmipaggregate,
    author = {Nicholls, Zebedee and Lewis, Jared and Makin, Melissa and Nattala, Usha and Zhang, Geordie Z. and Mutch, Simon J. and Tescari, Edoardo and Meinshausen, Malte},
    title = {Regionally aggregated, stitched and de-drifted CMIP-climate data, processed with netCDF-SCM {v2.0.0} (submitted)},
    journal = {Geoscience Data Journal},
    keywords = {climate, projections, model, CMIP, aggregate},
    year = {2020}
}

Licensing
=========

All data is shared under the same license terms as its source.
Again, to make life simpler, we provide simplified tools which provide a quick overview of your CMIP output's license terms (again see `our tools here <https://netcdf-scm.readthedocs.io/en/latest/usage/using-cmip-data.html>`_).

Retractions
===========

Occasionally CMIP data is retracted i.e. marked as erroneous.
We remove such data as quickly as possible.
However, if you have downloaded data, please check for retractions before using it and before publishing any conclusions based on it.
This can be done using the tools provided by netCDF-SCM (again see `<https://netcdf-scm.readthedocs.io/en/latest/usage/using-cmip-data.html>`_).
