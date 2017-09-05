Tool - Batch Correction [Dockerfile]
=======

Metadata
-----------

 * **@name**: Tool - Batch Correction [Dockerfile]
 * **@version**: 0.1
 * **@authors**: Nils Paulhe <nils.paulhe@inra.fr> (Only the docker part)
 * **@date creation**: 2017/09/04
 * **@main usage**: create a Docker environment / container for "Tool - Batch Correction"

 About
-----------
For all informations about the tool please refer to its [README file](README.txt). 
For further informations about W4M project and the people involved, Please refer to [W4M github](https://github.com/) and [W4M Docker Hub](https://hub.docker.com/r/workflow4metabolomics/). 
 
Configuration
-----------

### Requirement:
 * Docker Engine
 * Docker skills
 
### Warning:
 * These scripts are provided WITHOUT ANY WARRANTY. 
 * These scripts should be run by an administrator system (expert).

Services provided
-----------
Build a docker container for "Tool - Batch Correction" Galaxy Tool.
Provide a XML Galaxy wrapper: [batch_correction.docker.xml](batch_correction.docker.xml)

 
Technical description
-----------

### Create the docker container

``` bash
docker build -t workflow4metabolomics/tool-batch_correction:2.1.2 .
```

### Add the tool in Galaxy

Note: the files name and path are just examples. Adapt them to your own Galaxy configuration / practices.

If required, add in `config/job_conf.xml` file the minimal docker options:

``` xml
    <destinations default="docker_local">
        <destination id="local" runner="local"/>
        <destination id="docker_local" runner="local">
          <param id="docker_enabled">true</param>
          <param id="docker_sudo">false</param>
       </destination>
    </destinations>
```

For more options please refer to the [official documentation](https://galaxyproject.org/admin/tools/docker/).

Copy or create a symbolic link of [batch_correction.docker.xml](batch_correction.docker.xml) file into your `tools/docker` directory (feel free to create or change the target directory). 
Then add this XML resource in your `config/tool_conf.xml` file. For example:

``` xml
    <section id="docker_tools" name="Docker Tools">
      <tool file="docker/batch_correction.docker.xml"/>
    </section>
```

Notes
-----------
TODO

License (Dockerfile only!)
-----------
TODO