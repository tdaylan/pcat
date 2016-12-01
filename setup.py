from distutils.core import setup

setup( \
      name='pcat', \
      packages = ['pcat'], \
      version='1.0', \
      package_dir={'pcat': ''}, \
      packages=['pcat'], \
      description = 'A hierarchical, transdimensional MCMC sampler to explore the catalog space', \
      author = 'Tansu Daylan', \
      author_email = 'tansu.daylan@gmail.com', \
      download_url = 'https://github.com/tdaylan/pcat/tarball/0.1', \
      keywords = ['mcmc', 'bayesian', 'transdimensional', 'catalog', 'hierarchical'], \
     )

