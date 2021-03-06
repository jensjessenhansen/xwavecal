from setuptools import setup, find_packages


with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='xwavecal',
      author=['G. Mirek Brandt, Curtis McCully'],
      version='0.1.4',
      python_requires='>=3.6',
      packages=find_packages(),
      package_dir={'xwavecal': 'xwavecal'},
      package_data={'xwavecal': ['data/.*']},
      setup_requires=['pytest-runner'],
      install_requires=requirements,
      tests_require=['pytest>=3.5'],
      entry_points={'console_scripts': ['xwavecal_reduce_dir=xwavecal.main:run',
                                        'xwavecal_reduce=xwavecal.main:reduce_data']})
