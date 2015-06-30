from distutils.core import setup
from setuptools import setup

setup(
    name='umicount',
    version='1.0',
    packages=['umicount'],
    url='',
    license='GPL 2.0',
    author='mickael',
    author_email='mendez.mickael@gmail.com',
    description='A set a helpful function to count transcripts using UMI information ',

    entry_points = {
        'console_scripts': ['umicountFP=umicount.dedup_fingerprint:main'],
    }
)
