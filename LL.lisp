;; See Licence for legal information

(defparameter *cells* '())

(defclass std-cell ()
  ((POS
    :initarg :POS
    :accessor POS)
  (ATP
    :initarg :ATP
    :accessor ATP)
  (NA 
    :initarg :NA
    :accessor NA)
  (AA
    :initarg :AA
    :accessor AA)
  (FA
    :initarg :FA
    :accessor FA)
  (G
    :initarg :G
    :accessor G)
  (DNA
    :initarg :DNA
    :accessor DNA)))

(defun new-cell (POS ATP NA AA FA G DNA)
  (push (make-instance 'std-cell :POS POS :ATP ATP :NA NA :AA AA 
                                :FA FA :G G :DNA DNA) *cells*))

(defmacro fetch-value (accessor lst index)
  `(,accessor (nth ,index (reverse ,lst))))
