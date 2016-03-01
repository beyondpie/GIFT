import os
import sys

from pyplusplus import module_builder
mb = module_builder.module_builder_t(
    files=['hello_world.cpp'],
    gccxml_path='/home/zusongpeng/gccxml/bin')
mb.build_code_creator( module_name = 'libszu_sayhi_py')
mb.code_creator.user_defined_directories.append(os.path.abspath('.'))
mb.write_module(os.path.join(os.path.abspath('.'),'hello_world.cpp'))
