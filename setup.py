from skbuild import setup

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="Photochem",
    packages=['Photochem'],
    python_requires='>=3.7',
    version="0.0.1",
    license="MIT",
    install_requires=['numpy','scipy'],
    author='Nicholas Wogan',
    author_email = 'nicholaswogan@gmail.com',
    description = "Photochemical model of planet's atmospheres.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url = "https://github.com/Nicholaswogan/Photochem",
    include_package_data=True,
    cmake_args=['-DSKBUILD=ON',\
                '-DBUILD_PYTHON_PHOTOCHEM=ON',\
                '-DBUILD_STATIC_PHOTOCHEM=OFF']
)


