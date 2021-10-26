import setuptools

with open("readme.txt") as f:
    readme = f.read()

with open("LICENSE") as f:
    project_license = f.read()

with open("requirements.txt") as f:
    requirements = [line.strip() for line in f.readlines()]
    requirements_clean = [req.split('=')[0] for req in requirements]  # without '==version'

setuptools.setup(
    #
    name="AtlejgTools",
    package_dir = {'': 'lib/Python3'},
    packages=setuptools.find_packages('lib/Python3'),
    scripts=[ "lib/Python3/AtlejgTools/Scripts/run_pywake.py"],
    #
    install_requires=requirements,
    #install_requires=requirements_clean,
    #
    version="0.1.0",
    description="atle's stuff",
    long_description=readme,
    author="Equinor ASA",
    author_email="fg_sib-scout@equinor.com",
    url="https://github.com/atlejg/atlejg",
)
