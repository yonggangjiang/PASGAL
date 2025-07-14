# Root Makefile for MaxFlow

.PHONY: maxflow run_maxflow

maxflow:
	$(MAKE) -C src/MaxFlow

run_maxflow:
	@if [ -z "$(INPUT)" ]; then \
		echo "Please specify INPUT, e.g., make run_maxflow INPUT=your_input_file"; \
		exit 1; \
	fi; \
	./src/MaxFlow/maxflow $(INPUT)

# Usage:
#   make maxflow
#   make run_maxflow INPUT=your_input_file
