import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="galaxy-data-mines-mshubat",
    version="0.0.1",
    author="Matt Shubat",
    author_email="mattshubat@gmail.com",
    description="Galaxy Data Mine Application",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mshubat/galaxy_data_mines",
    packages=setuptools.find_packages(),
    package_dir={'galaxy_data_mines': 'galaxy_data_mines'},
    package_data={'galaxy_data_mines': ['data/*']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'gdmines=galaxy_data_mines.user_agent:main',
        ],
    },
)
