import setuptools

setuptools.setup(
    name="petace-duet",
    version="0.2.0",
    author="Tiktok PILab",
    author_email="",
    description="petace-duet",
    url="",
    install_requires=[
        "numpy"
    ],
    package_dir={"petace.duet": "./petace/duet/"},
    package_data={
        "petace.duet": ["*.so"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    packages=setuptools.find_packages(exclude=["petace.tests"]),
    python_requires="==3.9.*",
)
