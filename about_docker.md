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
For all informations about the tool please refer to its [README file](README.md). 
For further informations about Workflow4Metabolomics project and the people involved, please refer to [workflow4metabolomics.org](http://workflow4metabolomics.org/), [W4M github](https://github.com/workflow4metabolomics/) and [W4M Docker Hub](https://hub.docker.com/r/workflow4metabolomics/). 
 
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

License (Dockerfile only!)
-----------
The `Dockerfile` file is under the following license:
```
    Copyright (c) 2017 workflow4metabolomics.org / INRA

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
```

 * For the Galaxy Tool's license, please refer to its `README` file. 
 * For the Galaxy Wrapper's license, please refer to its `XML` file. 
