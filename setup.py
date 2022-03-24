from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()
setup(
    name='pyMultiOmics',
    version='1.0.4',
    author='Joe Wandy',
    author_email='joe.wandy@glasgow.ac.uk',
    description='A Python package for multi-omics data integration and analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/glasgowcompbio/pyMultiOmics',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    packages=find_packages(),
    package_data={
        'pyMultiOmics': [
            'data/*.*',
        ],
    },
    install_requires=['numpy', 'scipy', 'pandas', 'seaborn', 'scikit-learn', 'matplotlib', 'plotly',
                      'statsmodels', 'loguru', 'requests', 'neo4j-driver',
                      'tqdm', 'pillow', 'jupyterlab', 'bioservices', 'networkx', 'tzlocal', 'pals-pathway',
                      'mofapy2', 'mofax'],
)
