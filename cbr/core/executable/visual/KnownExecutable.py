from os import path
from PyQt5.QtCore import pyqtSignal, pyqtSlot, QThread
from PyQt5.QtWidgets import QDialog
import requests
import tempfile
from typing import Optional
import zipfile

from ..data import KnownExecutableLocation
from .Ui_KnownExecutable import Ui_KnownExecutable

class DownloadThread(QThread):

    progress_signal = pyqtSignal(object)

    def __init__(
        self,
        resource: str,
        location: str
    ) -> None:
        super().__init__()
        self.__resource = resource
        self.__location = location
        self.__error = None

    def error(self) -> Optional[Exception]:
        return self.__error

    def __download(self):

        with tempfile.TemporaryDirectory() as archive_dir:
            archive_file = path.join(archive_dir, "download.zip")
            r = requests.get(self.__resource, stream = True)
            size_str = r.headers.get('Content-Length')
            progress = 0

            if size_str is None:
                self.progress_signal.emit(None)
                size = None
            else:
                size = int(size_str)

            with open(archive_file, "wb") as archive:
                for chunk in r.iter_content(chunk_size = 1024):
                    if chunk:
                        progress += len(chunk)
                        archive.write(chunk)

                    if size is not None:
                        self.progress_signal.emit(progress / size)

            with zipfile.ZipFile(archive_file, 'r') as zip:
                zip.extractall(self.__location)

    def run(self):
        self.__error = None
        
        try:
            self.__download()
        except Exception as e:
            self.__error = e
            raise e

class KnownExecutableDialog(QDialog):

    def __init__(
        self,
        target : KnownExecutableLocation,
        components_folder : str
    ) -> None:
        super().__init__()
        self.__ui = Ui_KnownExecutable()
        self.__ui.setupUi(self)
        self.__download_thread = DownloadThread(
            target.location,
            components_folder
        )
        self.__download_thread.progress_signal.connect(self.__on_progress)
        self.__download_thread.finished.connect(self.__on_finished)

        self.__ui.infoTextEdit.setHtml(
            message.format(
                component = target.name,
                url = target.location,
                location = components_folder
            )
        )

        self.__ui.downloadButton.clicked.connect(self.__on_download_clicked)

    def error(self):
        return self.__download_thread.error()

    @pyqtSlot()
    def __on_finished(self):

        error = self.__download_thread.error()
        if error is None:
            self.accept()
        else:
            self.reject()

        self.__ui.downloadButton.setEnabled(True)

    @pyqtSlot(object)
    def __on_progress(self, value: Optional[float]):

        if value is None:
            return
        
        progress = self.__ui.progressBar
        top = progress.maximum()

        if top == 0:
            progress.setRange(0, 100)

        progress.setValue(int(value * 100))

    @pyqtSlot()
    def __on_download_clicked(self):
        progress = self.__ui.progressBar
        progress.setRange(0,0)
        self.__ui.downloadButton.setEnabled(False)
        self.__download_thread.start()

message ="""
The component <b>{component}</b> is needed to use this feature. You can click
the <b>Download</b> button which will automatically downlad the component from
<a href="{url}">{url}</a>. The component will be saved at <b>{location}</b>.
It can be deleted at anytime and this dialog will be displayed again if the
component is removed.
"""

