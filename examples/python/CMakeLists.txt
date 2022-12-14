foreach(example actor-create actor-daemon actor-join actor-kill actor-migrate actor-suspend actor-yield actor-lifetime
        app-masterworkers
        comm-wait comm-waitall comm-waitallfor comm-waitany comm-waitfor comm-failure comm-host2host comm-pingpong
        comm-ready comm-serialize comm-suspend comm-testany comm-throttling comm-waitallfor comm-waituntil
        exec-async exec-basic exec-dvfs exec-remote exec-ptask
        platform-profile platform-failures
        network-nonlinear clusters-multicpu io-degradation exec-cpu-nonlinear
        synchro-barrier synchro-mutex synchro-semaphore)
  set(tesh_files    ${tesh_files}   ${CMAKE_CURRENT_SOURCE_DIR}/${example}/${example}.tesh)
  set(examples_src  ${examples_src} ${CMAKE_CURRENT_SOURCE_DIR}/${example}/${example}.py)

  if(enable_python)
    ADD_TESH(python-${example}
                               --setenv srcdir=${CMAKE_CURRENT_SOURCE_DIR}/${example}
                               --setenv platfdir=${CMAKE_HOME_DIRECTORY}/examples/platforms
                               --setenv pythoncmd=${PYTHON_EXECUTABLE}
                               --setenv LD_LIBRARY_PATH=${TESH_LIBRARY_PATH}
                               --setenv PYTHONPATH=${CMAKE_BINARY_DIR}/lib
                               --cd ${CMAKE_CURRENT_SOURCE_DIR}/${example}
                               ${example}.tesh)
  endif()
endforeach()

set(examples_src  ${examples_src}                                                      PARENT_SCOPE)
set(tesh_files    ${tesh_files}   examples/python/actor-create/actor-create_d.xml
                                  examples/python/actor-lifetime/actor-lifetime_d.xml
                                  examples/python/app-masterworkers/app-masterworkers_d.xml
                                  examples/python/platform-failures/platform-failures_d.xml
                                  PARENT_SCOPE)
