module input_utils

  implicit none
  private
  public :: locate_block_list_file

contains

  subroutine locate_block_list_file(input_unit, block_list_file, has_file, input_filename)
    integer, intent(in) :: input_unit
    character(len=*), intent(out) :: block_list_file
    logical, intent(out) :: has_file
    character(len=*), intent(in), optional :: input_filename

    character(len=1024) :: actual_input_path
    character(len=1024) :: base_dir

    block_list_file = ''
    has_file = .false.

    if (present(input_filename) .and. len_trim(input_filename) > 0) then
      actual_input_path = trim(input_filename)
    else
      actual_input_path = ''
      call inquire_input_filename(input_unit, actual_input_path)
    end if

    if (len_trim(actual_input_path) == 0) return

    base_dir = parent_directory(actual_input_path)
    call find_block_list_reference(actual_input_path, base_dir, block_list_file, has_file)
    if (has_file .and. len_trim(block_list_file) > 0) then
      block_list_file = adjustl(trim(block_list_file))
    end if
  end subroutine locate_block_list_file

  subroutine inquire_input_filename(unit, filename)
    integer, intent(in) :: unit
    character(len=*), intent(out) :: filename
    character(len=1024) :: temp_name
    logical :: opened
    integer :: stat

    filename = ''
    inquire(unit=unit, name=temp_name, opened=opened, iostat=stat)
    if (opened .and. stat == 0) then
      filename = trim(temp_name)
    end if
  end subroutine inquire_input_filename

  subroutine find_block_list_reference(input_path, base_dir, block_list_file, has_file)
    character(len=*), intent(in) :: input_path, base_dir
    character(len=*), intent(out) :: block_list_file
    logical, intent(out) :: has_file

    integer :: scan_unit, stat
    character(len=1024) :: line, cleaned, block_line, candidate
    character(len=512) :: extracted
    integer :: pos_block, pos_eq

    has_file = .false.
    block_list_file = ''

    if (len_trim(input_path) == 0) return

    open(newunit=scan_unit, file=adjustl(input_path), status='old', action='read', iostat=stat)
    if (stat /= 0) return

    do
      read(scan_unit, '(A)', iostat=stat) line
      if (stat /= 0) exit
      cleaned = strip_comments(line)
      cleaned = adjustl(trim(cleaned))
      if (len_trim(cleaned) == 0) cycle
      pos_block = index(to_lower(cleaned), '&block_list')
      if (pos_block == 0) cycle

      block_line = cleaned(pos_block:)
      pos_eq = index(block_line, '=')
      if (pos_eq == 0) exit

      candidate = adjustl(block_line(pos_eq+1:))
      if (extract_quoted_token(candidate, extracted)) then
        block_list_file = make_full_path(base_dir, extracted)
        if (len_trim(block_list_file) > 0) then
          has_file = .true.
        end if
      end if
      exit
    end do

    close(scan_unit)
  end subroutine find_block_list_reference

  function strip_comments(text) result(clean)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: clean
    integer :: pos

    pos = index(text, '!')
    if (pos > 0) then
      clean = text(:pos-1)
    else
      clean = text
    end if
  end function strip_comments

  function to_lower(text) result(lowered)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lowered
    integer :: i, c

    lowered = text
    do i = 1, len(text)
      c = ichar(text(i:i))
      if (c >= ichar('A') .and. c <= ichar('Z')) then
        lowered(i:i) = achar(c + ichar('a') - ichar('A'))
      end if
    end do
  end function to_lower

  function extract_quoted_token(text, token) result(success)
    character(len=*), intent(in) :: text
    character(len=*), intent(out) :: token
    logical :: success
    integer :: start, finish
    character :: quote_char

    token = ''
    success = .false.

    start = index(text, "'")
    if (start /= 0) then
      quote_char = "'"
    else
      start = index(text, '"')
      if (start /= 0) then
        quote_char = '"'
      else
        return
      end if
    end if

    finish = index(text(start+1:), quote_char)
    if (finish == 0) return
    finish = start + finish

    token = trim(text(start+1:finish-1))
    if (len_trim(token) > 0) success = .true.
  end function extract_quoted_token

  function parent_directory(path) result(dir)
    character(len=*), intent(in) :: path
    character(len=len(path)) :: dir
    integer :: i, lenp

    dir = ''
    lenp = len_trim(path)
    if (lenp == 0) return

    do i = lenp, 1, -1
      if (path(i:i) == '/' .or. path(i:i) == '\\') then
        dir = path(:i)
        return
      end if
    end do
  end function parent_directory

  function make_full_path(base_dir, relative_path) result(full_path)
    character(len=*), intent(in) :: base_dir, relative_path
    character(len=1024) :: full_path
    character(len=1024) :: trimmed_base, trimmed_rel

    trimmed_rel = trim(relative_path)
    if (len_trim(trimmed_rel) == 0) then
      full_path = ''
      return
    end if

    if (is_absolute_path(trimmed_rel)) then
      full_path = trimmed_rel
      return
    end if

    trimmed_base = trim(base_dir)
    if (len_trim(trimmed_base) == 0) then
      full_path = trimmed_rel
      return
    end if

    if (trimmed_base(len_trim(trimmed_base):len_trim(trimmed_base)) /= '/') then
      trimmed_base = trimmed_base // '/'
    end if

    full_path = trimmed_base // trimmed_rel
  end function make_full_path

  function is_absolute_path(path) result(abs_path)
    character(len=*), intent(in) :: path
    logical :: abs_path

    abs_path = (len_trim(path) > 0) .and. (path(1:1) == '/')
  end function is_absolute_path

end module input_utils
