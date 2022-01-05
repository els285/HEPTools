import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="heptools",
    version="0.0.1",
    author="Ethan Simpson",
    author_email="ethansimpson285@gmail.com",
    description="Set of tools for high-energy physics analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ethansimpson285/HEPTools",
    project_urls={
        "Bug Tracker": "https://github.com/ethansimpson285/HEPTools/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)