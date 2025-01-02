# setup.py
import setuptools

setuptools.setup(
    name="amyg",
    version="0.1.0",             # adjust as needed
    author="Your Name",
    author_email="your_email@example.com",
    description="A pipeline for De Novo Genomic Annotation of Non-Model Organisms",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/cfarkas/amyg",   # your GitHub repo
    license="MIT",                           # matches your LICENSE
    py_modules=["amyg"],        # since you have amyg.py at top-level
    python_requires=">=3.7",
    # If you want a console script called "amyg" that calls amyg.py's main(),
    entry_points={
        "console_scripts": [
            "amyg=amyg:main",  
            # means: user can type 'amyg' in shell => runs amyg.py's main()
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
