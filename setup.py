import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "recognise",
    version = "0.0.1",
    author = "Nicholas Croucher",
    author_email = "nickjcroucher@imperial.ac.uk",
    description = ("Precise inference of bacterial recombinations"),
    license = "GPLv3",
    keywords = "bacterial genomics recombination bioinformatics",
    url = "https://github.com/nickjcroucher/recognise",
    packages=["recognise"],
    entry_points={
        "console_scripts": [
            "recognise = recognise.__main__:main"
        ]
    },
    test_suite="test",
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GPLv3 License",
    ],
)

