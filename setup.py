from setuptools import setup


setup(
    name="ugw",
    setup_requires=["setuptools_scm"],
    python_requires=">=3.6",
    install_requires=[
        "geopandas",
		"json",
		"owslib",
		"pandas",
		"pyproj",
		"requests",
    ],
    packages=["ugw"],
    package_dir={"ugw": "ugw"},
    zip_safe=False,
)
