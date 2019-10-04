from setuptools import setup

DISTUTILS_DEBUG = 1

config = dict(
    name='labber-script-extension',
    version='1.0',
    description='Extension for labber scripting tools',
    license="MIT",
    author='Antti Vepsalainen',
    author_email='avepsala@mit.edu',
    url='https://github.com/ittnas/labber-script-extension/',
    packages=['labber-script-extension'],
    # python_requires='>=3',
    install_requires=[]
)

setup(**config)
