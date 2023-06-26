from setuptools import setup, find_packages

setup(
    name='snakesuite',
    version='0.1.0',
    py_modules=['snakesuite'], 
    description='A toolbox of NGS workflows',
    url='https://github.com/niekwit/snakesuite',
    author='Niek Wit',
    author_email='nw416@cam.ac.uk',
    license='MIT',
    packages=find_packages(),
    install_requires=['seaborn','matplotlib','numpy','pandas',
                      'pyyaml','cutadapt', 'multiqc','Click'
                      ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License', 
        'Operating System :: POSIX :: Linux', 
        'Programming Language :: Python :: 3',
    ],
    entry_points={
        'console_scripts': [
            'snakesuite = snakesuite.scripts.snakesuite:cli',
        ],
    },
    include_package_data=True,
)