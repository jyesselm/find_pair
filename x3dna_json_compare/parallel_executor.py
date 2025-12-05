"""
Parallel execution utilities.

Provides wrappers for parallel processing of comparison tasks.
"""

import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Callable, List, Optional, TypeVar, Any
from functools import partial
import sys

T = TypeVar('T')
R = TypeVar('R')


class ParallelExecutor:
    """Wrapper for parallel execution with progress tracking."""
    
    def __init__(self, max_workers: Optional[int] = None, use_processes: bool = True):
        """
        Initialize parallel executor.
        
        Args:
            max_workers: Maximum number of workers (defaults to CPU count)
            use_processes: Use processes (True) or threads (False)
        """
        if max_workers is None:
            max_workers = max(1, multiprocessing.cpu_count() - 1)
        
        self.max_workers = max_workers
        self.use_processes = use_processes
        self.executor_class = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
    
    def execute_batch(self, tasks: List[Callable[[], R]],
                     progress_callback: Optional[Callable[[int, int], None]] = None) -> List[R]:
        """
        Execute tasks in parallel with progress tracking.
        
        Args:
            tasks: List of callable tasks (no arguments)
            progress_callback: Optional callback(count, total) for progress
            
        Returns:
            List of results in task order
        """
        results = [None] * len(tasks)
        
        with self.executor_class(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_index = {
                executor.submit(task): i
                for i, task in enumerate(tasks)
            }
            
            # Collect results as they complete
            completed = 0
            for future in as_completed(future_to_index):
                index = future_to_index[future]
                try:
                    results[index] = future.result()
                except Exception as e:
                    results[index] = e
                
                completed += 1
                if progress_callback:
                    progress_callback(completed, len(tasks))
        
        return results
    
    def map(self, func: Callable[[T], R], items: List[T],
           chunk_size: Optional[int] = None,
           progress_callback: Optional[Callable[[int, int], None]] = None) -> List[R]:
        """
        Parallel map with chunking for large datasets.
        
        Args:
            func: Function to apply to each item
            items: List of items to process
            chunk_size: Optional chunk size (for better progress tracking)
            progress_callback: Optional callback(count, total) for progress
            
        Returns:
            List of results in item order
        """
        if chunk_size is None:
            # Default to reasonable chunk size
            chunk_size = max(1, len(items) // (self.max_workers * 4))
        
        # Create tasks that process chunks
        def process_chunk(chunk: List[T]) -> List[R]:
            return [func(item) for item in chunk]
        
        chunks = [items[i:i+chunk_size] for i in range(0, len(items), chunk_size)]
        chunk_tasks = [partial(process_chunk, chunk) for chunk in chunks]
        
        chunk_results = self.execute_batch(
            chunk_tasks,
            progress_callback=lambda c, t: progress_callback(
                min(c * chunk_size, len(items)), len(items)
            ) if progress_callback else None
        )
        
        # Flatten results
        results = []
        for chunk_result in chunk_results:
            if isinstance(chunk_result, Exception):
                # Error occurred - re-raise or handle
                raise chunk_result
            results.extend(chunk_result)
        
        return results
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        return False


def print_progress(current: int, total: int):
    """Default progress printer."""
    percent = (current / total * 100) if total > 0 else 0
    sys.stdout.write(f"\rProgress: {current}/{total} ({percent:.1f}%)")
    sys.stdout.flush()
    if current == total:
        sys.stdout.write("\n")

