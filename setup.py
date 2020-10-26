from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="PALS-pathway",
    version="1.4.7",
    author="Joe Wandy",
    author_email="joe.wandy@glasgow.ac.uk",
    description="A Python tool to rank significantly-changing metabolite sets",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/glasgowcompbio/PALS",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': ['run_pals=pals.run_pals:main'],
    },
    scripts=['pals/run_pals.py'],
    python_requires='>=3.6',
    packages=find_packages(),
    package_data={
        'pals': [
            'data/*.json.zip',
            'data/reactome/*/*/*.json.zip'
        ],
    },
    install_requires=['numpy', 'scipy', 'pandas', 'seaborn', 'scikit-learn', 'matplotlib', 'statsmodels',
                      'loguru', 'requests', 'neo4j-driver', 'gseapy==0.9.16', 'streamlit', 'tqdm', 'pillow'],
)
