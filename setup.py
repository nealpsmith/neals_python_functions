from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))


setup(
    name='neals_python_functions',

    version='0.01',

    description='Repository of random functions that are useful',
    long_description='',

    url='',

    author='Neal Smith',
    author_email='nsmith19@mgh.harvard.edu',

    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'License :: OSI Approved :: MIT',

        'Programming Language :: Python :: 3',
    ],

    install_requires=["numpy"],

    keywords='',

    packages=find_packages(),

    entry_points={
        'console_scripts': [
            'sample=sample:main',
        ],
    },
)