import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="neop",
    version="0.1.0",
    author="Ananthan Sadagopan",
    author_email="ananthans1000@gmail.com",
    description="Python Package for Predicting Neoantigens and Obtaining Amino Acid Context from MAFs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/as1000/neop",
    project_urls={
        "Bug Tracker": "https://github.com/as1000/neop/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    py_modules=['main'],
    python_requires=">=3.7",
)
