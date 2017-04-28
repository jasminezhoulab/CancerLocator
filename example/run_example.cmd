java -jar ..\CancerLocator.jar example_config > run.log
FC run_result reference_result.dos >NUL && Echo Test run passed || Echo Test run failed
