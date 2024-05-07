from os import path
import os
from PyQt5.QtCore import QDir, QUrl
from PyQt5.QtWidgets import QWidget
from tempfile import TemporaryDirectory
from typing import Any, Callable, Dict, List, TypeVar

from .context.visual.RunQWidgetContext import RunQWidgetContext
from .executable import *

T = TypeVar('T')

MAIN_USER_MANUAL = QUrl('https://github.com/TUM-CBR/pymol-plugins/wiki')

AppCloseHandler = Callable[[], Any]

class Context(object):

    def __init__(
        self,
        user_manual: Optional[QUrl] = None
    ):
        self.__store : Dict[str, Any] = {}
        self.__widgets : List[RunQWidgetContext] = []
        self.__directories : List[TemporaryDirectory[Any]] = []
        self.__application_directory = path.join(
            QDir.homePath(),
            ".cbrtools"
        )
        self.__user_manual = user_manual if user_manual is not None else MAIN_USER_MANUAL
        self.__child_contexts: List['Context'] = []
        self.__app_close_handlers: List[AppCloseHandler] = []

    def __ensure_directory_exists(self, directory: str) -> str:
        if path.exists(directory):
            return directory
        
        os.mkdir(directory)
        return directory

    def __get_application_directory(self) -> str:
        return self.__ensure_directory_exists(self.__application_directory)

    def __components_directory(self):
        return self.__ensure_directory_exists(
            path.join(self.__get_application_directory(), "components")
        )
    
    def get_config_directory(self, name: str):

        dir = path.join(
            self.__get_application_directory(),
            name
        )

        return self.__ensure_directory_exists(dir)

    def create_temporary_directory(self) -> str:
        directory = TemporaryDirectory()
        self.__directories.append(directory)
        return directory.name
    
    def __cleanup__(self):
        for directory in self.__directories:
            directory.cleanup()

        self.__directories.clear()

        for ctx in self.__child_contexts:
            ctx.__cleanup__()

        self.__child_contexts.clear()

    def __del__(self):
        self.__cleanup__()

    def get_executable(
        self,
        widget: QWidget,
        executable: KnownExecutables
    ) -> Executable:
        
        return find_executable(
            widget,
            executable,
            self.__components_directory()
        )
                

    def create_or_load(self, name : str, load : Callable[[], T]) -> T:

        if name in self.__store:
            return self.__store[name]

        result = self.__store[name] = load()

        return result

    def on_app_close(self, handler: AppCloseHandler):
        self.__app_close_handlers.append(handler)

    def __on_close(self, window: RunQWidgetContext):

        for handler in self.__app_close_handlers:

            try:
                handler()
            except Exception:
                # Handler threw excepotion, we need to
                # continue execution
                pass

        self.__app_close_handlers = []

        try:
            self.__widgets.remove(window)
        except ValueError:
            # This should not happen, but not critical if it does
            pass

    def create_child_context(
        self,
        new_user_manual: Optional[QUrl] = None
    ) -> 'Context':
        
        ctx = Context(
            user_manual = new_user_manual if new_user_manual is not None else self.__user_manual
        )
        self.__child_contexts.append(ctx)

        return ctx
    
    def close_app(self):
        for widet in self.__widgets:
            widet.close()

    def run_widget(
        self,
        run_widget : Callable[['Context'], QWidget]
    ) -> RunQWidgetContext:

        widget = run_widget(self)
        window = RunQWidgetContext(widget, self.__on_close, self.__user_manual)

        self.__widgets.append(window)

        return window