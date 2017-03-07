java -jar ../CancerLocator.jar example_config > run.log
cmp -s run_result reference_result && echo "Test run passed" || echo "Test run failed"
