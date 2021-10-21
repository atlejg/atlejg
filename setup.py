import setuptools

with open("readme.txt") as f:
    readme = f.read()

with open("LICENSE") as f:
    project_license = f.read()

with open("requirements.txt") as f:
    requirements = [line.strip() for line in f.readlines()]

setuptools.setup(
    name="atlejg",
    version="0.1.0",
    description="atle's stuff",
    long_description=readme,
    author="Equinor ASA",
    author_email="fg_sib-scout@equinor.com",
    url="https://github.com/atlejg/atlejg",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    scripts=[ "lib/Python3/AtlejgTools/Scripts/run_pywake.py"],
)
