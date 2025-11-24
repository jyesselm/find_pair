# Makefile for x3dna project
# Provides convenient targets for common build and development tasks

.PHONY: help release debug clean format test install all check rebuild org-release org-debug org-relwithdebinfo org-clean org-clean-all

# Default target
.DEFAULT_GOAL := help

# Detect OS
UNAME_S := $(shell uname -s)

# Detect number of CPU cores
ifeq ($(UNAME_S),Darwin)
    # macOS
    NUM_CPUS := $(shell sysctl -n hw.ncpu)
else
    # Linux
    NUM_CPUS := $(shell nproc)
endif

# Directories
BUILD_DIR := build
PROJECT_ROOT := $(shell pwd)

# Build system (prefer Ninja if available)
HAS_NINJA := $(shell command -v ninja 2> /dev/null)
ifdef HAS_NINJA
    GENERATOR := Ninja
    BUILD_CMD := ninja -j $(NUM_CPUS)
else
    GENERATOR := Unix Makefiles
    BUILD_CMD := $(MAKE) -j$(NUM_CPUS)
endif

# Colors for output (if terminal supports it)
ifneq ($(TERM),)
    COLOR_RESET := \033[0m
    COLOR_BOLD := \033[1m
    COLOR_GREEN := \033[32m
    COLOR_YELLOW := \033[33m
    COLOR_BLUE := \033[34m
endif

help: ## Show this help message
	@echo "$(COLOR_BOLD)Available targets:$(COLOR_RESET)"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  $(COLOR_GREEN)%-15s$(COLOR_RESET) %s\n", $$1, $$2}'
	@echo ""
	@echo "$(COLOR_BOLD)Examples:$(COLOR_RESET)"
	@echo "  make release     # Build optimized release version"
	@echo "  make debug       # Build with debug symbols"
	@echo "  make test        # Run all tests"
	@echo "  make format      # Format all C++ source files"
	@echo "  make clean       # Remove build artifacts"

release: ## Build optimized release version (-O3, -DNDEBUG)
	@echo "$(COLOR_BOLD)Building Release...$(COLOR_RESET)"
	@$(MAKE) _build BUILD_TYPE=Release

debug: ## Build with debug symbols (-g, -O0)
	@echo "$(COLOR_BOLD)Building Debug...$(COLOR_RESET)"
	@$(MAKE) _build BUILD_TYPE=Debug

relwithdebinfo: ## Build optimized with debug symbols (-O2, -g)
	@echo "$(COLOR_BOLD)Building RelWithDebInfo...$(COLOR_RESET)"
	@$(MAKE) _build BUILD_TYPE=RelWithDebInfo

_build: ## Internal build target (don't call directly)
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && \
	if [ ! -f CMakeCache.txt ] || \
	   [ "$(PROJECT_ROOT)/CMakeLists.txt" -nt CMakeCache.txt ] || \
	   [ "$$(grep -i 'CMAKE_BUILD_TYPE:STRING' CMakeCache.txt 2>/dev/null | cut -d'=' -f2)" != "$(BUILD_TYPE)" ] || \
	   [ "$$(grep -i 'CMAKE_GENERATOR:INTERNAL' CMakeCache.txt 2>/dev/null | cut -d'=' -f2 | grep -i $(GENERATOR) || echo '')" = "" ]; then \
		echo "$(COLOR_YELLOW)Configuring CMake with $(GENERATOR) generator...$(COLOR_RESET)"; \
		cmake "$(PROJECT_ROOT)" \
			-G "$(GENERATOR)" \
			-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
			-DCMAKE_EXPORT_COMPILE_COMMANDS=ON; \
	fi
	@cd $(BUILD_DIR) && $(BUILD_CMD)
	@echo "$(COLOR_GREEN)Build complete!$(COLOR_RESET)"

clean: ## Remove build artifacts (keeps CMakeCache)
	@echo "$(COLOR_BOLD)Cleaning build artifacts...$(COLOR_RESET)"
	@if [ -d $(BUILD_DIR) ]; then \
		cd $(BUILD_DIR) && $(BUILD_CMD) clean 2>/dev/null || true; \
	fi
	@echo "$(COLOR_GREEN)Clean complete!$(COLOR_RESET)"

clean-all: ## Remove entire build directory
	@echo "$(COLOR_BOLD)Removing build directory...$(COLOR_RESET)"
	@rm -rf $(BUILD_DIR)
	@echo "$(COLOR_GREEN)Clean complete!$(COLOR_RESET)"

format: ## Format all C++ source files using clang-format
	@echo "$(COLOR_BOLD)Formatting C++ source files...$(COLOR_RESET)"
	@find . -type f \( -name "*.cpp" -o -name "*.hpp" -o -name "*.h" -o -name "*.cc" -o -name "*.cxx" \) \
		! -path "./$(BUILD_DIR)/*" \
		! -path "./.git/*" \
		! -path "./org/build/*" \
		! -path "./Testing/*" \
		-exec clang-format -i {} + 2>/dev/null || \
		$(MAKE) format-script
	@echo "$(COLOR_GREEN)Formatting complete!$(COLOR_RESET)"

format-script: ## Format using format_code.sh script
	@if [ -f scripts/format_code.sh ]; then \
		bash scripts/format_code.sh; \
	else \
		echo "$(COLOR_YELLOW)Warning: clang-format not found and format_code.sh not available$(COLOR_RESET)"; \
	fi

format-check: ## Check code formatting without modifying files
	@echo "$(COLOR_BOLD)Checking code formatting...$(COLOR_RESET)"
	@find . -type f \( -name "*.cpp" -o -name "*.hpp" -o -name "*.h" -o -name "*.cc" -o -name "*.cxx" \) \
		! -path "./$(BUILD_DIR)/*" \
		! -path "./.git/*" \
		! -path "./org/build/*" \
		! -path "./Testing/*" \
		-exec sh -c 'clang-format "$$1" | diff -u "$$1" - || exit 1' _ {} \; || \
		(echo "$(COLOR_YELLOW)Formatting issues found. Run 'make format' to fix.$(COLOR_RESET)" && exit 1)
	@echo "$(COLOR_GREEN)Formatting check passed!$(COLOR_RESET)"

test: ## Run all tests
	@echo "$(COLOR_BOLD)Running tests...$(COLOR_RESET)"
	@if [ ! -d $(BUILD_DIR) ]; then \
		echo "$(COLOR_YELLOW)Build directory not found. Building in Debug mode first...$(COLOR_RESET)"; \
		$(MAKE) debug; \
	fi
	@cd $(BUILD_DIR) && ctest --output-on-failure -j$(NUM_CPUS)
	@echo "$(COLOR_GREEN)Tests complete!$(COLOR_RESET)"

test-verbose: ## Run tests with verbose output
	@echo "$(COLOR_BOLD)Running tests (verbose)...$(COLOR_RESET)"
	@if [ ! -d $(BUILD_DIR) ]; then \
		echo "$(COLOR_YELLOW)Build directory not found. Building in Debug mode first...$(COLOR_RESET)"; \
		$(MAKE) debug; \
	fi
	@cd $(BUILD_DIR) && ctest --output-on-failure --verbose -j$(NUM_CPUS)

rebuild: clean-all release ## Clean and rebuild from scratch in release mode

rebuild-debug: clean-all debug ## Clean and rebuild from scratch in debug mode

install: ## Install the project (if installation is configured)
	@if [ ! -d $(BUILD_DIR) ]; then \
		echo "$(COLOR_YELLOW)Build directory not found. Building in Release mode first...$(COLOR_RESET)"; \
		$(MAKE) release; \
	fi
	@cd $(BUILD_DIR) && cmake --install . || echo "$(COLOR_YELLOW)Installation not configured$(COLOR_RESET)"

check: format-check test ## Run format check and tests

all: release ## Build release version (default build target)

compile_commands: ## Generate compile_commands.json for IDE support
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && \
	cmake "$(PROJECT_ROOT)" \
		-G "$(GENERATOR)" \
		-DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_EXPORT_COMPILE_COMMANDS=ON
	@if [ -f $(BUILD_DIR)/compile_commands.json ]; then \
		ln -sf $(BUILD_DIR)/compile_commands.json . 2>/dev/null || true; \
		echo "$(COLOR_GREEN)compile_commands.json generated!$(COLOR_RESET)"; \
	else \
		echo "$(COLOR_YELLOW)Warning: compile_commands.json not found$(COLOR_RESET)"; \
	fi

run-tool: ## Run generate_modern_json tool (usage: make run-tool ARGS="<pdb_file> <json_file> [--legacy]")
	@if [ ! -f $(BUILD_DIR)/generate_modern_json ]; then \
		echo "$(COLOR_YELLOW)Tool not found. Building in Release mode first...$(COLOR_RESET)"; \
		$(MAKE) release; \
	fi
	@$(BUILD_DIR)/generate_modern_json $(ARGS)

org-release: ## Build org code in Release mode
	@echo "$(COLOR_BOLD)Building org code (Release)...$(COLOR_RESET)"
	@mkdir -p org/build
	@cd org/build && \
	if [ -f CMakeCache.txt ]; then \
		CACHED_GEN=$$(grep -i 'CMAKE_GENERATOR:INTERNAL' CMakeCache.txt 2>/dev/null | cut -d'=' -f2 | tr -d ' '); \
		if [ -n "$$CACHED_GEN" ] && ! echo "$$CACHED_GEN" | grep -qi "$(GENERATOR)"; then \
			echo "$(COLOR_YELLOW)Generator mismatch. Cleaning org build directory...$(COLOR_RESET)"; \
			cd .. && rm -rf build && mkdir -p build && cd build; \
		fi; \
	fi && \
	if [ ! -f CMakeCache.txt ] || \
	   [ "../CMakeLists.txt" -nt CMakeCache.txt ] || \
	   [ "$$(grep -i 'CMAKE_BUILD_TYPE:STRING' CMakeCache.txt 2>/dev/null | cut -d'=' -f2)" != "Release" ]; then \
		echo "$(COLOR_YELLOW)Configuring CMake for org code...$(COLOR_RESET)"; \
		cmake .. -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=Release; \
	fi
	@cd org/build && $(BUILD_CMD)
	@echo "$(COLOR_GREEN)Org code build complete!$(COLOR_RESET)"
	@echo "  Executables: org/build/bin/find_pair_analyze, find_pair_original, analyze_original"

org-debug: ## Build org code in Debug mode
	@echo "$(COLOR_BOLD)Building org code (Debug)...$(COLOR_RESET)"
	@mkdir -p org/build
	@cd org/build && \
	if [ -f CMakeCache.txt ]; then \
		CACHED_GEN=$$(grep -i 'CMAKE_GENERATOR:INTERNAL' CMakeCache.txt 2>/dev/null | cut -d'=' -f2 | tr -d ' '); \
		if [ -n "$$CACHED_GEN" ] && ! echo "$$CACHED_GEN" | grep -qi "$(GENERATOR)"; then \
			echo "$(COLOR_YELLOW)Generator mismatch. Cleaning org build directory...$(COLOR_RESET)"; \
			cd .. && rm -rf build && mkdir -p build && cd build; \
		fi; \
	fi && \
	if [ ! -f CMakeCache.txt ] || \
	   [ "../CMakeLists.txt" -nt CMakeCache.txt ] || \
	   [ "$$(grep -i 'CMAKE_BUILD_TYPE:STRING' CMakeCache.txt 2>/dev/null | cut -d'=' -f2)" != "Debug" ]; then \
		echo "$(COLOR_YELLOW)Configuring CMake for org code...$(COLOR_RESET)"; \
		cmake .. -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=Debug; \
	fi
	@cd org/build && $(BUILD_CMD)
	@echo "$(COLOR_GREEN)Org code build complete!$(COLOR_RESET)"
	@echo "  Executables: org/build/bin/find_pair_analyze, find_pair_original, analyze_original"

org-relwithdebinfo: ## Build org code in RelWithDebInfo mode
	@echo "$(COLOR_BOLD)Building org code (RelWithDebInfo)...$(COLOR_RESET)"
	@mkdir -p org/build
	@cd org/build && \
	if [ -f CMakeCache.txt ]; then \
		CACHED_GEN=$$(grep -i 'CMAKE_GENERATOR:INTERNAL' CMakeCache.txt 2>/dev/null | cut -d'=' -f2 | tr -d ' '); \
		if [ -n "$$CACHED_GEN" ] && ! echo "$$CACHED_GEN" | grep -qi "$(GENERATOR)"; then \
			echo "$(COLOR_YELLOW)Generator mismatch. Cleaning org build directory...$(COLOR_RESET)"; \
			cd .. && rm -rf build && mkdir -p build && cd build; \
		fi; \
	fi && \
	if [ ! -f CMakeCache.txt ] || \
	   [ "../CMakeLists.txt" -nt CMakeCache.txt ] || \
	   [ "$$(grep -i 'CMAKE_BUILD_TYPE:STRING' CMakeCache.txt 2>/dev/null | cut -d'=' -f2)" != "RelWithDebInfo" ]; then \
		echo "$(COLOR_YELLOW)Configuring CMake for org code...$(COLOR_RESET)"; \
		cmake .. -G "$(GENERATOR)" -DCMAKE_BUILD_TYPE=RelWithDebInfo; \
	fi
	@cd org/build && $(BUILD_CMD)
	@echo "$(COLOR_GREEN)Org code build complete!$(COLOR_RESET)"
	@echo "  Executables: org/build/bin/find_pair_analyze, find_pair_original, analyze_original"

org-clean: ## Clean org build artifacts (keeps CMakeCache)
	@echo "$(COLOR_BOLD)Cleaning org build artifacts...$(COLOR_RESET)"
	@if [ -d org/build ]; then \
		if [ -f org/build/build.ninja ]; then \
			cd org/build && ninja clean 2>/dev/null || true; \
		elif [ -f org/build/Makefile ]; then \
			cd org/build && $(MAKE) clean 2>/dev/null || true; \
		fi; \
	fi
	@echo "$(COLOR_GREEN)Org clean complete!$(COLOR_RESET)"

org-clean-all: ## Remove entire org build directory
	@echo "$(COLOR_BOLD)Removing org build directory...$(COLOR_RESET)"
	@rm -rf org/build
	@echo "$(COLOR_GREEN)Org clean complete!$(COLOR_RESET)"

info: ## Show build system information
	@echo "$(COLOR_BOLD)Build System Information:$(COLOR_RESET)"
	@echo "  OS:           $(UNAME_S)"
	@echo "  CPU cores:    $(NUM_CPUS)"
	@echo "  Generator:    $(GENERATOR)"
	@echo "  Build dir:    $(BUILD_DIR)"
	@echo "  Build cmd:    $(BUILD_CMD)"
	@echo ""
	@echo "$(COLOR_BOLD)CMake version:$(COLOR_RESET)"
	@cmake --version | head -1 || echo "  CMake not found"
	@echo ""
	@echo "$(COLOR_BOLD)clang-format:$(COLOR_RESET)"
	@clang-format --version 2>/dev/null | head -1 || echo "  clang-format not found"
	@if [ -f $(BUILD_DIR)/CMakeCache.txt ]; then \
		echo ""; \
		echo "$(COLOR_BOLD)Current build configuration:$(COLOR_RESET)"; \
		grep -i "CMAKE_BUILD_TYPE:STRING" $(BUILD_DIR)/CMakeCache.txt 2>/dev/null | cut -d'=' -f2 || echo "  Not configured"; \
	fi

