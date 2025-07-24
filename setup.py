from setuptools import setup, find_packages

setup(
    name='virus_variant_calling_pipeline',
    version='1.0.0',
    author='Rajindra Napit',
    description='Pipeline for dengue virus variant calling',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline',
    packages=find_packages(exclude=['conda-recipe']),
    install_requires=['pandas>=2.3.0','matplotlib','numpy','pyyaml'],
    entry_points={
        'console_scripts': [
            'virus-pipeline=run_pipeline:main'
        ]
    },
    python_requires='>=3.8',
)
