#!/usr/bin/env python
# -*- coding: utf-8 -*-


from distutils.core import setup

exec(open('info.py').read())


def readme():
    with open('README.md') as f:
        return f.read()

python_names = ['run_pipeline.py', 'cfis_download_images.py', 'cfis_field_select.py', 'find_mask_flag.py',
                  'scp_CFIS_cc.py', 'cfis_check_weights.py', 'cfis_get_coord_exposures.py',
                  'cfis_create_exposures.py', 'create_image_links.py', 'cfis_write_tileobj_as_exposures.py',
                  'cfis_select_tileobj_expPSF.py', 'test_mexp_CFIS_MOBJ.py']
sh_names     = ['remove_hdu0.sh']

scripts_py     = ['python/{}'.format(fn) for fn in python_names]
scripts_sh     = ['sh/{}'.format(fn) for fn in sh_names]

setup(name             = 'scripts',
      author           = __author__,
      packages         = ['cfis', 'generic'],
      package_dir      = {'cfis': 'python/lib/cfis', 'generic': 'python/lib/generic'},
      version          = __version__,
      author_email     = __email__,
      license          = 'GNU GPLv3',
      description      = 'Scripts and modules that can be used independently of the pipeline',
      platforms        = ['posix', 'mac os'],
      long_description = readme(),
      classifiers      = [
                          'Programming Language :: Python',
                          'Natural Language :: English',
                         ],
      scripts          = scripts_py + scripts_sh,
    )
 
