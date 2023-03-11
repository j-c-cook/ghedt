from setuptools import setup, find_packages

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='ghedt',
      install_requires=['pygfunction>=2.1',
                        'pandas>=1.3.2',
                        'openpyxl>=3.0.8',
                        'opencv-python==4.5.4.58'],
      url='https://github.com/j-c-cook/ghedt',
      download_url='https://github.com/j-c-cook/ghedt/archive/v0.3.1.tar.gz',
      long_description=long_description,
      long_description_content_type='text/markdown',
      version='0.3.1',
      packages=find_packages(),
      include_package_data=True,
      author='Jack C. Cook',
      author_email='jack-c-cook@protonmail.com',
      description='A ground heat exchanger design tool with the advanced and '
                  'unmatched capability of automatic borehole field selection '
                  'based on drilling geometric land constraints.')
