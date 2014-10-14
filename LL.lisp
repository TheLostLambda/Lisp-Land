(defvar *cell-id* 0)

(defclass std-cell ()
  ((ID
    :initform (incf *cell-id*))
  (POS
    :initarg :POS
    :initform (cons 100 100))
  (ATP
    :initarg :ATP
    :initform 1000)
  (NA 
    :initarg :NA
    :initform 1000)
  (AA
    :initarg :AA
    :initform 1200000)
  (FA
    :initarg :FA
    :initform 5000)
  (G
    :initarg :G
    :initform 100)
  (DNA
    :initarg :DNA
    :initform (error "Genome must be specified"))))
