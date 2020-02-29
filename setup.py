import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="OOpyPST-marinlauber", # Replace with your own username
    version="0.0.1",
    author="Marin Lauber",
    author_email="M.Lauber@soton.ac.uk",
    description="OOP 2D pseudo-spectral code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marinlauber/OOpyPST",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
)
