<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>main</name>
    <title>MPF - Multiprocessing framework</title>
    <filename>main</filename>
    <docanchor file="main">Description</docanchor>
    <docanchor file="main">Execution</docanchor>
    <docanchor file="main">Installation</docanchor>
    <docanchor file="main">notes</docanchor>
    <docanchor file="main">Prerequisites</docanchor>
    <docanchor file="main">code</docanchor>
  </compound>
  <compound kind="file">
    <name>__init__.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>____init_____8py</filename>
    <namespace>mpfg</namespace>
  </compound>
  <compound kind="file">
    <name>mp_args.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__args_8py</filename>
    <class kind="class">mpfg::mp_args::Args</class>
    <class kind="class">mpfg::mp_args::Args::ParsingError</class>
    <class kind="class">mpfg::mp_args::Args::ValidationError</class>
    <namespace>mpfg::mp_args</namespace>
  </compound>
  <compound kind="file">
    <name>mp_calc.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__calc_8py</filename>
    <class kind="class">mpfg::mp_calc::Master</class>
    <class kind="class">mpfg::mp_calc::Worker</class>
    <namespace>mpfg::mp_calc</namespace>
  </compound>
  <compound kind="file">
    <name>mp_calc_MPI.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__calc__MPI_8py</filename>
    <class kind="class">mpfg::mp_calc_MPI::MasterMPI</class>
    <class kind="class">mpfg::mp_calc_MPI::WorkerMPI</class>
    <namespace>mpfg::mp_calc_MPI</namespace>
    <member kind="variable">
      <type>int</type>
      <name>TAG_WORK</name>
      <anchorfile>namespacempfg_1_1mp__calc__MPI.html</anchorfile>
      <anchor>aae95ddedfd7b208aedd58e9d17ac5a80</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>TAG_RESULT</name>
      <anchorfile>namespacempfg_1_1mp__calc__MPI.html</anchorfile>
      <anchor>a9c3a397e2fad5bec5b7e9930579b57aa</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>TAG_TERMINATE</name>
      <anchorfile>namespacempfg_1_1mp__calc__MPI.html</anchorfile>
      <anchor>a296e7f48452321439fef14707cd65791</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mp_calc_SMP.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__calc__SMP_8py</filename>
    <class kind="class">mpfg::mp_calc_SMP::MasterSMP</class>
    <class kind="class">mpfg::mp_calc_SMP::WorkerSMP</class>
    <namespace>mpfg::mp_calc_SMP</namespace>
  </compound>
  <compound kind="file">
    <name>mp_data.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__data_8py</filename>
    <class kind="class">mpfg::mp_data::Dataset</class>
    <namespace>mpfg::mp_data</namespace>
  </compound>
  <compound kind="file">
    <name>mp_helper.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__helper_8py</filename>
    <class kind="class">mpfg::mp_helper::Helper</class>
    <namespace>mpfg::mp_helper</namespace>
  </compound>
  <compound kind="file">
    <name>mp_job.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__job_8py</filename>
    <class kind="class">mpfg::mp_job::JobProcessor</class>
    <class kind="class">mpfg::mp_job::Job</class>
    <class kind="class">mpfg::mp_job::FileJob</class>
    <class kind="class">mpfg::mp_job::JobResult</class>
    <namespace>mpfg::mp_job</namespace>
  </compound>
  <compound kind="file">
    <name>mp_MPI.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__MPI_8py</filename>
    <namespace>mpfg::mp_MPI</namespace>
    <member kind="variable">
      <type></type>
      <name>comm</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>aac59683b65e2a8113476ea99b23b38ff</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>rank</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>a6e43ae43e56ea61488a95a6b92ef7cd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>args</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>a6b1e9c92d030d6de7960ad5ff167275e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>master</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>aa7ebc9ab8415854c5e11e282f0192f05</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>worker</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>a9d6ca8202732490d76c570df94bcf0f2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mp_SMP.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__SMP_8py</filename>
    <namespace>mpfg::mp_SMP</namespace>
    <member kind="variable">
      <type>tuple</type>
      <name>args</name>
      <anchorfile>namespacempfg_1_1mp__SMP.html</anchorfile>
      <anchor>aa924536eb6d9894971428c9c44015c93</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>master</name>
      <anchorfile>namespacempfg_1_1mp__SMP.html</anchorfile>
      <anchor>a47d4622a002df11435333b53484b1263</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mp_version.py</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>mp__version_8py</filename>
    <namespace>mpfg::mp_version</namespace>
    <namespace>mpfg::version</namespace>
    <member kind="variable">
      <type>string</type>
      <name>__version__</name>
      <anchorfile>namespacempfg_1_1mp__version.html</anchorfile>
      <anchor>a9774249d56277f1904f579e735581841</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Exception</name>
    <filename>classException.html</filename>
  </compound>
  <compound kind="class">
    <name>object</name>
    <filename>classobject.html</filename>
  </compound>
  <compound kind="namespace">
    <name>mpfg</name>
    <filename>namespacempfg.html</filename>
    <namespace>mpfg::mp_args</namespace>
    <namespace>mpfg::mp_calc</namespace>
    <namespace>mpfg::mp_calc_MPI</namespace>
    <namespace>mpfg::mp_calc_SMP</namespace>
    <namespace>mpfg::mp_data</namespace>
    <namespace>mpfg::mp_helper</namespace>
    <namespace>mpfg::mp_job</namespace>
    <namespace>mpfg::mp_MPI</namespace>
    <namespace>mpfg::mp_SMP</namespace>
    <namespace>mpfg::mp_version</namespace>
    <namespace>mpfg::version</namespace>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_args</name>
    <filename>namespacempfg_1_1mp__args.html</filename>
    <class kind="class">mpfg::mp_args::Args</class>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_args::Args</name>
    <filename>classmpfg_1_1mp__args_1_1Args.html</filename>
    <base>object</base>
    <class kind="class">mpfg::mp_args::Args::ParsingError</class>
    <class kind="class">mpfg::mp_args::Args::ValidationError</class>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a0c727e3cc2a66a79311159e3a07538b5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>helper</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a4eac2d056b5b2a9764e20e115dd19716</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>args</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a0650902e55fcc49b77f19d08ad427a3d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>options</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a52a7779294be2d473208c57e6e331962</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>default_config_dir</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a8078c09325305ca7f780e494cdcfa606</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>default_config_filename</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a24ac1a6d3922ec7de373bd64bd726b6c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>args</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a0650902e55fcc49b77f19d08ad427a3d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>options</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a52a7779294be2d473208c57e6e331962</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>default_config_dir</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a8078c09325305ca7f780e494cdcfa606</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>default_config_filename</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a24ac1a6d3922ec7de373bd64bd726b6c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>parse_command_line</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>aa54dae0f52aff4188b236962506fd661</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>parse_input_arguments</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>afe571fb318ca26672e6cd9f795d5c881</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>validate_arguments</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>ab5c4de3143161186829fbd2fcf5cd9be</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>add_option</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a8ca0e4aff91faa12394a516cbd692502</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>print_usage</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a2a00ff6d36aa80ad7e4e73143d7e6382</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>get_prog_name</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a827e42171f0edd076ffb2af8fac0018c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type></type>
      <name>args</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a3ec549c4880ae33c1eb143616fe79a1d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>string</type>
      <name>err_msg</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a7fb540a0c0972fe27cd5f82afe358b48</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type></type>
      <name>validated</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>aa8ece017244c917117d143fea3abb588</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>tuple</type>
      <name>config_filepath</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args.html</anchorfile>
      <anchor>a921dcfff039f2458c347e5aafddb0b4f</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_args::Args::ParsingError</name>
    <filename>classmpfg_1_1mp__args_1_1Args_1_1ParsingError.html</filename>
    <base>Exception</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args_1_1ParsingError.html</anchorfile>
      <anchor>aa80838f82268e6538776b9407ccbb542</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>__str__</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args_1_1ParsingError.html</anchorfile>
      <anchor>aa8a181e5ff858ccff62204a657e90042</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_args::Args::ValidationError</name>
    <filename>classmpfg_1_1mp__args_1_1Args_1_1ValidationError.html</filename>
    <base>Exception</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args_1_1ValidationError.html</anchorfile>
      <anchor>a8aeb1b310b311d8ae97502a9e6c3a01f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>__str__</name>
      <anchorfile>classmpfg_1_1mp__args_1_1Args_1_1ValidationError.html</anchorfile>
      <anchor>a2a4ad5c2fb3acfb5646b4060e92c478d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_calc</name>
    <filename>namespacempfg_1_1mp__calc.html</filename>
    <class kind="class">mpfg::mp_calc::Master</class>
    <class kind="class">mpfg::mp_calc::Worker</class>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_calc::Master</name>
    <filename>classmpfg_1_1mp__calc_1_1Master.html</filename>
    <base>object</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a640e1cea426f594558b54e2631786c4b</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>__str__</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a343a6d521e0a157fffb3a0f5032993b4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>helper</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>aa5bbcd6990529c186d60b72e9e06cdfe</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>arg_options</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>aca2df551b2af6ad135f95f3d761966c8</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>config</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a6198ee998ef81d3de54984635a33d35a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>default_config_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a86dc2b39b7dd205097060321b6abdf6c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>default_config_filename</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a112dc5dd78b106ae32781abcce95b17e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>config_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>aaf7f846f9d3536c6ea57686f4d7aed9f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>config_filename</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a46c20d143814dc90e5ae6e51f6c8d660</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>logger</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>ae210ded69762f1dec2ff669f98062803</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>name</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a3bd21497d723b97c89dc3cec0c2d79e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>base_input_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a4c7b7cd03b49ae85dfdf4a018ae713ba</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>base_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>aa49814afa9fe5759b4b0dd8848b9032f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>abe144b72993b0ff5313c8fdda6b4e27b</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>log_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a6dc878cf2d19e7057c6584652df18d6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>plot_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a10baee0a8d8e98fd7cabbd2272ee36ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>map_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a855451a067398db81a6eb13824316c93</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>stat_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>ab7c05ca581d8852b075994fb7afc5245</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>result_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a7cc521d0198da63b5fa982b974b0ebdb</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>error_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a60e8cffd1fe7d43f9d62b537dca94af5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>job_processor</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>ada3ca016fae3a3026b81d6e7eb4652da</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>job_list</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a7c08a4d9d8446b91f8e287fff173e755</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>shared_dico</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a742a3fa776f19ac449b0610790e42c1a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run_duration</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>ac43840688aad1fedfb2fa22895494c26</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>name</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a3bd21497d723b97c89dc3cec0c2d79e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>job_processor</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>ada3ca016fae3a3026b81d6e7eb4652da</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>logger</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>ae210ded69762f1dec2ff669f98062803</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>logging_enabled</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a5c0ba3c4dd0229042d8867f854631c9f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a8d46c92ac98570881fb22e0de09b01d6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>shutdown</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a9e848a1dd1c08ef3d57a4bb09c4d2df4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>process_job_result</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a24be0396ba974ad3cd520c975ff098fa</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>get_default_config_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a33da4f87ebdc0ac857627bdd993db185</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>get_default_config_filename</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>aa5fa0d877e864e76419cca7e8908ff38</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>create_logger</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a7b2e90e5b74d2294e2efd050a872b3d9</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>create_run_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>aa9be0e8988ebc140d4787f3f03ad6114</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>create_job_processor</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a1a2cf6d8eb1a53f92edca1bcf6b98f22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type></type>
      <name>logger</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Master.html</anchorfile>
      <anchor>a1e47a47c19e53a0e69ef80461b622b40</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_calc::Worker</name>
    <filename>classmpfg_1_1mp__calc_1_1Worker.html</filename>
    <base>object</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>abe8cadc06beb989cd8f0c96cc991d7b2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>__str__</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>abb04edf838a8a8adc5986dfa97dd184e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>helper</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a2ff6b98c4fe32d6cf116faa9f3094160</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>arg_options</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a08b5eeffa54d97db1eb7cb1aed7228e6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>config</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a921397bb185079bd937a1f706e4ae34c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>default_config_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>ad5bbe44db27290210d54b2b05bb1b667</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>logger</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a6e87bad87a22b3ff5be6143ec26877dc</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>name</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a8b7778b4768b296223bc4ed562cf8f35</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>shared_dico</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>aeb5e36f245ae0d0ed0821c0138ad3f34</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>job_processor</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a6e0a332e16857bbae4ed8e9eed13874c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run_duration</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>aee467517b2ef155e9ac5f705ad1cc171</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>base_input_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a38dbd8290aff0e69c769c7711c181dfb</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>base_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a31b61c6efc10a70b41fecf4e48f36695</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a6d1e574a2bf65a4ef492b2767579618a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>log_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>ab66bf86ef49411ae4cb7005e70fb7e45</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>plot_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a612a1b132141e1b636c53bc8de1661c1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>map_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a8660eb2e51172306e9468b96875706d5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>stat_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a0883fdeccb04d25f70066eb9dadb59db</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>result_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a4010f2c8c17bd442793a2cfaaeae2291</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>error_output_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>ac44004ec7a1e0abdfb4d61c5e8509e5f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>logger</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a6e87bad87a22b3ff5be6143ec26877dc</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>name</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a8b7778b4768b296223bc4ed562cf8f35</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>shared_dico</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>aeb5e36f245ae0d0ed0821c0138ad3f34</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>get_default_config_dir</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>aea85079e189586e81f2a1da1d8a138ac</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>get_default_config_filename</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>aa149e586a4e678a6b72b7bb74a5685e8</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>logging_enabled</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a2f1af7edd9af86c7d8af3f3ecded28dc</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>preprocess_job</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>aa5f8d9677cdc1239e298c4bea709820d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>process_job</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a62fb4745285dc341123883af55eedb6c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>postprocess_job</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>acb236396db2abcf9bc3a07558afbd161</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>create_logger</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>a8af44b7ec05bb89e3b2d57f018829818</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type></type>
      <name>logger</name>
      <anchorfile>classmpfg_1_1mp__calc_1_1Worker.html</anchorfile>
      <anchor>abb8dfde36c21611774dca4e39daa48f8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_calc_MPI</name>
    <filename>namespacempfg_1_1mp__calc__MPI.html</filename>
    <class kind="class">mpfg::mp_calc_MPI::MasterMPI</class>
    <class kind="class">mpfg::mp_calc_MPI::WorkerMPI</class>
    <member kind="variable">
      <type>int</type>
      <name>TAG_WORK</name>
      <anchorfile>namespacempfg_1_1mp__calc__MPI.html</anchorfile>
      <anchor>aae95ddedfd7b208aedd58e9d17ac5a80</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>TAG_RESULT</name>
      <anchorfile>namespacempfg_1_1mp__calc__MPI.html</anchorfile>
      <anchor>a9c3a397e2fad5bec5b7e9930579b57aa</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>TAG_TERMINATE</name>
      <anchorfile>namespacempfg_1_1mp__calc__MPI.html</anchorfile>
      <anchor>a296e7f48452321439fef14707cd65791</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_calc_MPI::MasterMPI</name>
    <filename>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</filename>
    <base>mpfg::mp_calc::Master</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</anchorfile>
      <anchor>a4581221df4170fbcdab418dbaddd721e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>comm</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</anchorfile>
      <anchor>a0c4a1ce007f72dc0a52b455e0215bc50</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>rank</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</anchorfile>
      <anchor>a534f9fc11cd442a53308fbdc48b996a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>nb_workers</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</anchorfile>
      <anchor>aaa97894526a29d2a943ab39ec69c50f3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>nb_workers</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</anchorfile>
      <anchor>aaa97894526a29d2a943ab39ec69c50f3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>input_queue</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</anchorfile>
      <anchor>acf1ff7147dcfa6e63df76a37109ba64d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>submit_jobs</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</anchorfile>
      <anchor>a2725bdbd65bd908074f9f12db7d981f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</anchorfile>
      <anchor>abba30ba7f190fa137e57955c5b27a60e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>process_jobs</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1MasterMPI.html</anchorfile>
      <anchor>ab43e509e2fdb8913123f21067b9f4172</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_calc_MPI::WorkerMPI</name>
    <filename>classmpfg_1_1mp__calc__MPI_1_1WorkerMPI.html</filename>
    <base>mpfg::mp_calc::Worker</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1WorkerMPI.html</anchorfile>
      <anchor>a296d5ece6d6b4e7f465738577897672b</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>comm</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1WorkerMPI.html</anchorfile>
      <anchor>aee632f4e18ba84937f768485132e8f9d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>rank</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1WorkerMPI.html</anchorfile>
      <anchor>a1af0a1ca11df42709c85c5655bacc4b3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1WorkerMPI.html</anchorfile>
      <anchor>a081d251c69b56535c8677b1431ce61a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type></type>
      <name>shared_dico</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1WorkerMPI.html</anchorfile>
      <anchor>a3c4d3cc868425b879cf540713225b4ed</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type></type>
      <name>logger</name>
      <anchorfile>classmpfg_1_1mp__calc__MPI_1_1WorkerMPI.html</anchorfile>
      <anchor>abcd272f8b3f36a68c4d920cb7d38cc7d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_calc_SMP</name>
    <filename>namespacempfg_1_1mp__calc__SMP.html</filename>
    <class kind="class">mpfg::mp_calc_SMP::MasterSMP</class>
    <class kind="class">mpfg::mp_calc_SMP::WorkerSMP</class>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_calc_SMP::MasterSMP</name>
    <filename>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</filename>
    <base>mpfg::mp_calc::Master</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>acd219350b15e6cdad4ccd2a8948c1006</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>workers</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>a4885fa4f27f6285feb52012ee50c20c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>nb_workers</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>aa7424fbe9c085252af312e8c89c1c6ee</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>input_queue</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>a25d0930c097efb294b1d844a5a81e4f8</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>output_queue</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>af87f9308625e402bdd4097f8b9458916</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>start_event</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>afc9e3a1eacea9d8ab8a1f4655934b456</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>workers</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>a4885fa4f27f6285feb52012ee50c20c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>a0d764152591f9cda89813f6941894a80</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>create_workers</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>a85014b1a7d5c30e44eec598bf4516c72</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>create_worker</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>a30a12f39ff62ba822b875c32e85913e5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>submit_jobs</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>a12d8fc9c6336805c5e830109329d81ed</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run_workers</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>a228d643360c4dbd0643402a44c85c9d1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>process_jobs</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>ad747898324051320a08457bf2fd36faa</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>all_jobs_processed</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1MasterSMP.html</anchorfile>
      <anchor>a1f4dae6ae18da403fb7d51561670ac23</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_calc_SMP::WorkerSMP</name>
    <filename>classmpfg_1_1mp__calc__SMP_1_1WorkerSMP.html</filename>
    <base>mpfg::mp_calc::Worker</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1WorkerSMP.html</anchorfile>
      <anchor>a63da1e6610ce57aa40c0cbbff11a6735</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>run</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1WorkerSMP.html</anchorfile>
      <anchor>a87c37d266405cb6f7cc0ab7b043d442f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type></type>
      <name>name</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1WorkerSMP.html</anchorfile>
      <anchor>a1c80e82d59bd38b9631e417a411a6511</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type></type>
      <name>shared_dico</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1WorkerSMP.html</anchorfile>
      <anchor>afcd28da9adf65f2ac3590de07a41bafd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type></type>
      <name>logger</name>
      <anchorfile>classmpfg_1_1mp__calc__SMP_1_1WorkerSMP.html</anchorfile>
      <anchor>a39e976b152f970c4f45cf9887aa349a2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_data</name>
    <filename>namespacempfg_1_1mp__data.html</filename>
    <class kind="class">mpfg::mp_data::Dataset</class>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_data::Dataset</name>
    <filename>classmpfg_1_1mp__data_1_1Dataset.html</filename>
    <base>object</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>a6843e4506638cd333eee06d1ebe87472</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>__str__</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>a0dfed3bb8f58bf43cdf1969377b64bd2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>helper</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>aa77db1a085c959e0e49f41f3c2e58a84</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>type</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>ab56ba04561c89419b87abb5190d50636</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>name</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>ad6e9e89ab77b27191b852ee583897b54</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>base_dir</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>a47f1158bad1968c3061945c5c2bb61f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>base_dir</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>a47f1158bad1968c3061945c5c2bb61f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>query</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>acb78ede5e36e46f239e80164f82e4dfb</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>locate_files</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>abcb2166829ab129c5126f22c77e30a7b</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>match_file</name>
      <anchorfile>classmpfg_1_1mp__data_1_1Dataset.html</anchorfile>
      <anchor>a34b24691ad35a51754e1f8026cfcb5b3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_helper</name>
    <filename>namespacempfg_1_1mp__helper.html</filename>
    <class kind="class">mpfg::mp_helper::Helper</class>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_helper::Helper</name>
    <filename>classmpfg_1_1mp__helper_1_1Helper.html</filename>
    <base>object</base>
    <member kind="function">
      <type>def</type>
      <name>file_exists</name>
      <anchorfile>classmpfg_1_1mp__helper_1_1Helper.html</anchorfile>
      <anchor>a6058e212696e60b12eb0f4e019da1885</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>get_main_basename</name>
      <anchorfile>classmpfg_1_1mp__helper_1_1Helper.html</anchorfile>
      <anchor>a991e3a28fc4acd55196a03a61bb4d1bc</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>read_config</name>
      <anchorfile>classmpfg_1_1mp__helper_1_1Helper.html</anchorfile>
      <anchor>aa447d475b55020a59f6ffe5ae9c5cc15</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>dump_config</name>
      <anchorfile>classmpfg_1_1mp__helper_1_1Helper.html</anchorfile>
      <anchor>a68e47789f5c1404537ed0ca2955d2596</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>is_logging</name>
      <anchorfile>classmpfg_1_1mp__helper_1_1Helper.html</anchorfile>
      <anchor>a8ca8344b3538bd169361e589e14fa419</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>make_dir</name>
      <anchorfile>classmpfg_1_1mp__helper_1_1Helper.html</anchorfile>
      <anchor>a996f925557302697bac46c38c87d9d3c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>print_error</name>
      <anchorfile>classmpfg_1_1mp__helper_1_1Helper.html</anchorfile>
      <anchor>a04b4ebc45468a653986b38bf086eceab</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>print_warning</name>
      <anchorfile>classmpfg_1_1mp__helper_1_1Helper.html</anchorfile>
      <anchor>a5df6522441b58a6a3cb69fd682f80543</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>print_info</name>
      <anchorfile>classmpfg_1_1mp__helper_1_1Helper.html</anchorfile>
      <anchor>a3057b1167608692042d85a76d457504e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_job</name>
    <filename>namespacempfg_1_1mp__job.html</filename>
    <class kind="class">mpfg::mp_job::JobProcessor</class>
    <class kind="class">mpfg::mp_job::Job</class>
    <class kind="class">mpfg::mp_job::FileJob</class>
    <class kind="class">mpfg::mp_job::JobResult</class>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_job::JobProcessor</name>
    <filename>classmpfg_1_1mp__job_1_1JobProcessor.html</filename>
    <base>object</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>aa42599fe3d0cb7380d5dc3f82366897f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>helper</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>a309776ca73c09be88d1f45515df4a1c0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>job_list</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>aa4499d8249cef1e14c30500e6e24b632</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>job_result_list</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>a557109652f59a2cb17c9ae53ab4c473a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>job_count</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>a26548d398d306aecbb84ab4647b9453d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>preprocess_job</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>ab09fcfd1722b883cc5786f9e1444d11c</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>process_job</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>ac82e281585aea3994c19f62a2c9c5072</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>postprocess_job</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>a09a93db2e165da3ab33dbd03db1aecb5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>process_job_result</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>afcd0d991640ce796a21decd051c2af09</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>create_jobs</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>a762d5ea58fee88b9b04fb36229435ee2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>create_dataset</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>a0e8d863863bf9ad8aecec9560b8a4615</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>record_job_result</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>a941392337bff164eca16f72124f6f3f3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>all_jobs_processed</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobProcessor.html</anchorfile>
      <anchor>a4531fb9bedd66330525b17739bc960e1</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_job::Job</name>
    <filename>classmpfg_1_1mp__job_1_1Job.html</filename>
    <base>object</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__job_1_1Job.html</anchorfile>
      <anchor>a79cdf9d86b50a3c8c25f0e63d3d19beb</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>__str__</name>
      <anchorfile>classmpfg_1_1mp__job_1_1Job.html</anchorfile>
      <anchor>acaae44f3573346d1d6ee96fdcc234981</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>name</name>
      <anchorfile>classmpfg_1_1mp__job_1_1Job.html</anchorfile>
      <anchor>a4d0f0b23d0f8bf05a9cb2425d5a8fb27</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>dataset</name>
      <anchorfile>classmpfg_1_1mp__job_1_1Job.html</anchorfile>
      <anchor>a57f799307835e73a5adad5d5ec747c31</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>name</name>
      <anchorfile>classmpfg_1_1mp__job_1_1Job.html</anchorfile>
      <anchor>a4d0f0b23d0f8bf05a9cb2425d5a8fb27</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_job::FileJob</name>
    <filename>classmpfg_1_1mp__job_1_1FileJob.html</filename>
    <base>mpfg::mp_job::Job</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__job_1_1FileJob.html</anchorfile>
      <anchor>ab1ee503332394cc27fb8a9dedf5f0bcc</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>filepath</name>
      <anchorfile>classmpfg_1_1mp__job_1_1FileJob.html</anchorfile>
      <anchor>a5605b6b6189ce3b84204d1326c792d88</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>filename</name>
      <anchorfile>classmpfg_1_1mp__job_1_1FileJob.html</anchorfile>
      <anchor>aa6dcf9d55de7b341ea04c11572434356</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>directory</name>
      <anchorfile>classmpfg_1_1mp__job_1_1FileJob.html</anchorfile>
      <anchor>ae4687914f7c1ea83c0164b8e5cea3dec</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>filepath</name>
      <anchorfile>classmpfg_1_1mp__job_1_1FileJob.html</anchorfile>
      <anchor>a5605b6b6189ce3b84204d1326c792d88</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>mpfg::mp_job::JobResult</name>
    <filename>classmpfg_1_1mp__job_1_1JobResult.html</filename>
    <base>object</base>
    <member kind="function">
      <type>def</type>
      <name>__init__</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobResult.html</anchorfile>
      <anchor>a2f3c9fe8ac9f1a24d8f40ae180ff6266</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>job</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobResult.html</anchorfile>
      <anchor>a36fdf0987d00942b467e9acc8bad6c92</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>result</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobResult.html</anchorfile>
      <anchor>afa8c7b96541b44bfc5958e687aafe190</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>def</type>
      <name>result</name>
      <anchorfile>classmpfg_1_1mp__job_1_1JobResult.html</anchorfile>
      <anchor>afa8c7b96541b44bfc5958e687aafe190</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_MPI</name>
    <filename>namespacempfg_1_1mp__MPI.html</filename>
    <member kind="variable">
      <type></type>
      <name>comm</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>aac59683b65e2a8113476ea99b23b38ff</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>rank</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>a6e43ae43e56ea61488a95a6b92ef7cd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>args</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>a6b1e9c92d030d6de7960ad5ff167275e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>master</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>aa7ebc9ab8415854c5e11e282f0192f05</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>worker</name>
      <anchorfile>namespacempfg_1_1mp__MPI.html</anchorfile>
      <anchor>a9d6ca8202732490d76c570df94bcf0f2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_SMP</name>
    <filename>namespacempfg_1_1mp__SMP.html</filename>
    <member kind="variable">
      <type>tuple</type>
      <name>args</name>
      <anchorfile>namespacempfg_1_1mp__SMP.html</anchorfile>
      <anchor>aa924536eb6d9894971428c9c44015c93</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>tuple</type>
      <name>master</name>
      <anchorfile>namespacempfg_1_1mp__SMP.html</anchorfile>
      <anchor>a47d4622a002df11435333b53484b1263</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::mp_version</name>
    <filename>namespacempfg_1_1mp__version.html</filename>
    <member kind="variable">
      <type>string</type>
      <name>__version__</name>
      <anchorfile>namespacempfg_1_1mp__version.html</anchorfile>
      <anchor>a9774249d56277f1904f579e735581841</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mpfg::version</name>
    <filename>namespacempfg_1_1version.html</filename>
  </compound>
  <compound kind="dir">
    <name>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</name>
    <path>/home/marc/cea-epfl/trunk/mpfg_package/mpfg/</path>
    <filename>dir_a1c9c5cd095ac3d5ae5ae551082fe7dc.html</filename>
    <file>__init__.py</file>
    <file>mp_args.py</file>
    <file>mp_calc.py</file>
    <file>mp_calc_MPI.py</file>
    <file>mp_calc_SMP.py</file>
    <file>mp_data.py</file>
    <file>mp_helper.py</file>
    <file>mp_job.py</file>
    <file>mp_MPI.py</file>
    <file>mp_SMP.py</file>
    <file>mp_version.py</file>
  </compound>
</tagfile>
