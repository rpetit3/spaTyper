import setuptools
import glob

with open("README", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="spaTyper",
    version="0.2.2",
    scripts=glob.glob('main/*'),
    author="Mitchell Sullivan; Jose F. Sanchez-Herrero",
    author_email="jfbioinformatics@gmail.com",
    description="Typing of Staphylococcus aureus protein A from assemblies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JFsanchezherrero/spa_typing",
    packages=setuptools.find_packages(),
    license='GPLv3',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
