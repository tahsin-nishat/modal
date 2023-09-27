# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

# Parse the version from the module.
with open('sysid/__init__.py') as f:
    for line in f:
        if line.find("__version__") >= 0:
            version = line.split("=")[1].strip()
            version = version.strip('"')
            version = version.strip("'")
            break

with open('README.rst', 'r', encoding="utf-8") as fin:
    long_description = fin.read()


setup(
    name='sysid',
    version=version,
    author='Tahsin Nishat',
    author_email='nishat@arizona.edu',
    description='System identification of a simply supported beam',
    license='ISC',
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url='https://github.com/tahsin-nishat/sysid',
    packages=find_packages(exclude=["test"]),
    classifiers=[
        "License :: OSI Approved :: ISC License (ISCL)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        ],
    install_requires=["numpy", "scipy", "matplotlib"],
    python_requires=">=3.6",
)
