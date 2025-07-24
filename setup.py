from setuptools import setup, find_packages

setup(
    name='virus_variant_calling_pipeline',
    version='1.0.0',
    author='Rajindra Napit',
    author_email='rajindra@d4@gmail.com',
    description='Pipeline for dengue virus variant calling and annotation',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline',
    packages=find_packages(exclude=['conda-recipe', 'tests']),
    include_package_data=True,
    install_requires=[
        'pandas>=2.3.0',
        'numpy>=2.3.1',
        'matplotlib>=3.10.3',
        'pyyaml>=6.0.2'
    ],
    entry_points={
        'console_scripts': [
            'virus-pipeline=run_pipeline:main'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
    ],
    python_requires='>=3.8',
)
