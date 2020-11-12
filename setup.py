import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='Praxis',
    version="0.2",
    url='https://github.com/davidlevybooth/Praxis',
    license='GPL-3',
    author='David Levy-Booth, Parker Lloyd',
    author_email='dlevyboo@mail.ubc.ca',
    description='Praxis - Prokaryotic Activity and Expression Infomatics System',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['Praxis'],
    package_data={'': [
            "Praxis/*",
                       ]},
    data_files=[(".", ["README.md", "LICENSE.txt"])],
    include_package_data=True,
    install_requires= [
    ],
    # install via conda: click, pandas, pyyaml, snakemake
    entry_points={
          'console_scripts': [
              'Praxis = Praxis.Praxis:main'
          ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)