from .designabc import DesignABC


class Design(DesignABC):
    def precompile(self):
        self.uncompile()
        with self.logger.timeit("DEBUG", "running blast"):
            self._blast()

    def postcompile(self, post_process_kwargs: dict = None):
        self.post_process_graphs(post_process_kwargs)
