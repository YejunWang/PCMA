from setuptools import setup, find_packages

setup(
    name='PCMA',
    version='0.1.0',
    packages=find_packages(),
    package_data={'pcma': ['rcode/heatmap.R']},
    install_requires=[
        'numpy', 'pandas', 'scikit-learn', 'scipy', 'statsmodels', 'joblib',
        'kaleido', 'jinja2', 'weasyprint', 'pdf2image', 'plotly', 'PyPDF2'
    ],
    author='yukki & Eathon',
    description='A package for PCA mediation analysis',
    long_description=open('readme.md').read(),
    long_description_content_type='text/markdown',
    url='NA',
)
